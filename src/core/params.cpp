#include "core/params.hpp"

#include <cmath>
#include <filesystem>
#include <limits>
#include "third_party/toml.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

// 一次 list-directed "记录读取"的机械复刻：READ(10,*) 无变量 = 丢弃一行；
// 带变量 READ = 取一行切 token。行数不足/文件缺失抛错。
class LineCursor {
public:
    explicit LineCursor(const std::string& path) : path_(path), in_(path) {
        if (!in_) throw std::runtime_error("Params: cannot open " + path);
    }

    void skip() {              // READ(10,*) 无变量：丢弃一行
        (void)next_line();
    }

    std::vector<std::string> tokens() {   // 带变量 READ：取一行切 token
        return split(next_line());
    }

private:
    std::string next_line() {
        std::string line;
        if (!std::getline(in_, line))
            throw std::runtime_error("Params: unexpected EOF in " + path_ +
                                     " at line " + std::to_string(lineno_ + 1));
        ++lineno_;
        return line;
    }

    static std::vector<std::string> split(const std::string& line) {
        // list-directed：逗号与空白都是分隔符
        std::vector<std::string> out;
        std::string cur;
        for (char c : line) {
            if (c == ',' || c == ' ' || c == '\t' || c == '\r') {
                if (!cur.empty()) { out.push_back(cur); cur.clear(); }
            } else {
                cur.push_back(c);
            }
        }
        if (!cur.empty()) out.push_back(cur);
        return out;
    }

    std::string path_;
    std::ifstream in_;
    int lineno_ = 0;
};

[[noreturn]] void bad_token(const std::string& file, const std::string& tok) {
    throw std::runtime_error("Params: bad numeric token '" + tok + "' in " + file);
}

double to_double(const std::string& file, const std::string& tok) {
    std::string t = tok;
    for (char& c : t)
        if (c == 'D' || c == 'd') c = 'E';
    try {
        size_t pos;
        double v = std::stod(t, &pos);
        if (pos != t.size()) bad_token(file, tok);
        return v;
    } catch (const std::invalid_argument&) { bad_token(file, tok); }
      catch (const std::out_of_range&)    { bad_token(file, tok); }
}

int to_int(const std::string& file, const std::string& tok) {
    // int 变量遇小数点/指数报错（§8：防隐式类型错位）
    if (tok.find_first_of(".eEdD") != std::string::npos)
        throw std::runtime_error("Params: integer expected, got '" + tok +
                                 "' in " + file);
    try {
        size_t pos;
        int v = std::stoi(tok, &pos);
        if (pos != tok.size()) bad_token(file, tok);
        return v;
    } catch (const std::invalid_argument&) { bad_token(file, tok); }
      catch (const std::out_of_range&)    { bad_token(file, tok); }
}

void expect(const std::string& file, const std::vector<std::string>& tok, size_t n) {
    if (tok.size() != n)
        throw std::runtime_error("Params: expected " + std::to_string(n) +
                                 " values, got " + std::to_string(tok.size()) +
                                 " in " + file + " (line: " + [&] {
                                     std::string s;
                                     for (auto& t : tok) s += t + " ";
                                     return s;
                                 }() + ")");
}

std::string join_path(const std::string& dir, const std::string& name) {
    return dir.empty() ? name : (dir.back() == '/' ? dir + name : dir + "/" + name);
}


// Fold ASCII configuration identifiers without changing string or numeric values.
std::string lowercase(std::string_view text) {
    std::string result(text);
    for (char& c : result)
        if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
    return result;
}

toml::table normalize_keys(const toml::table& source, const std::string& path = {}) {
    toml::table result;
    for (const auto& [key, value] : source) {
        const auto name = lowercase(key.str());
        const auto full_name = path.empty() ? name : path + "." + name;
        if (result.contains(name))
            throw std::runtime_error("duplicate TOML name (case-insensitive): " + full_name);
        if (const auto* child = value.as_table())
            result.insert(name, normalize_keys(*child, full_name));
        else
            result.insert(name, value);
    }
    return result;
}

Params load_toml(const std::string& path) {
    const auto table = normalize_keys(toml::parse_file(path));
    Params p{};
    auto section = [&](const char* name, std::initializer_list<const char*> keys) {
        auto t = table[name].as_table();
        if (!t) throw std::runtime_error(std::string("missing TOML table: ") + name);
        for (auto&& [key, value] : *t) {
            bool known = false;
            for (auto k : keys) known |= key.str() == lowercase(k);
            if (!known) throw std::runtime_error("unknown TOML key: " + std::string(name) + "." + std::string(key.str()));
        }
    };
    for (auto&& [key, value] : table)
        if (key != "mesh" && key != "parm" && key != "control" && key != "output")
            throw std::runtime_error("unknown TOML table: " + std::string(key.str()));
    auto integer = [&](const char* group, const char* key) {
        const auto n = table[group][lowercase(key)].value<int64_t>();
        if (!table[group][lowercase(key)].is_integer() || !n || *n < std::numeric_limits<int>::min() || *n > std::numeric_limits<int>::max())
            throw std::runtime_error(std::string("integer required: ") + group + "." + key);
        return static_cast<int>(*n);
    };
    auto real = [&](const char* group, const char* key) {
        const auto n = table[group][lowercase(key)].value<double>();
        if (!n || !std::isfinite(*n)) throw std::runtime_error(std::string("finite number required: ") + group + "." + key);
        return *n;
    };
    section("mesh", {"NX", "NY", "NZ", "INDEX_disturb_Z", "NXS", "NYS", "NZS", "NX2", "NY2", "NZ2", "INDEX_Serial"});
    p.mesh.NX = integer("mesh", "NX");
    p.mesh.NY = integer("mesh", "NY");
    p.mesh.NZ = integer("mesh", "NZ");
    p.mesh.INDEX_disturb_Z = integer("mesh", "INDEX_disturb_Z");
    p.mesh.NXS = integer("mesh", "NXS");
    p.mesh.NYS = integer("mesh", "NYS");
    p.mesh.NZS = integer("mesh", "NZS");
    p.mesh.NX2 = integer("mesh", "NX2");
    p.mesh.NY2 = integer("mesh", "NY2");
    p.mesh.NZ2 = integer("mesh", "NZ2");
    p.mesh.INDEX_Serial = integer("mesh", "INDEX_Serial");
    section("parm", {"INDEX_GRAD_PRE", "INDEX_UNSTEADY_BC", "XL", "ZL", "THETA", "Re", "G_pre", "a_pre0", "b_pre0", "a_V", "b_V"});
    p.parm.INDEX_GRAD_PRE = integer("parm", "INDEX_GRAD_PRE");
    p.parm.INDEX_UNSTEADY_BC = integer("parm", "INDEX_UNSTEADY_BC");
    p.parm.XL = real("parm", "XL");
    p.parm.ZL = real("parm", "ZL");
    p.parm.THETA = real("parm", "THETA");
    p.parm.Re = real("parm", "Re");
    p.parm.G_pre = real("parm", "G_pre");
    p.parm.a_pre0 = real("parm", "a_pre0");
    p.parm.b_pre0 = real("parm", "b_pre0");
    {
        auto a = table["parm"]["a_v"].as_array();
        if (!a || a->size() != 3) throw std::runtime_error("parm.a_V requires three numbers");
        for (size_t i = 0; i < 3; ++i) {
            auto v = (*a)[i].value<double>();
            if (!v || !std::isfinite(*v)) throw std::runtime_error("invalid parm.a_V element");
            p.parm.a_V[i] = *v;
        }
    }
    {
        auto a = table["parm"]["b_v"].as_array();
        if (!a || a->size() != 3) throw std::runtime_error("parm.b_V requires three numbers");
        for (size_t i = 0; i < 3; ++i) {
            auto v = (*a)[i].value<double>();
            if (!v || !std::isfinite(*v)) throw std::runtime_error("invalid parm.b_V element");
            p.parm.b_V[i] = *v;
        }
    }
    section("control", {"INDEX_CONT", "NSTEPS", "NDUMFRQ", "NPRNFRQ", "NCOUNT_disturb", "NDISTURBFRQ", "Nstat", "EPS_init", "DT", "EPS_disturb", "V_Microblow0"});
    p.ctrl.INDEX_CONT = integer("control", "INDEX_CONT");
    p.ctrl.NSTEPS = integer("control", "NSTEPS");
    p.ctrl.NDUMFRQ = integer("control", "NDUMFRQ");
    p.ctrl.NPRNFRQ = integer("control", "NPRNFRQ");
    p.ctrl.NCOUNT_disturb = integer("control", "NCOUNT_disturb");
    p.ctrl.NDISTURBFRQ = integer("control", "NDISTURBFRQ");
    p.ctrl.Nstat = integer("control", "Nstat");
    p.ctrl.EPS_init = real("control", "EPS_init");
    p.ctrl.DT = real("control", "DT");
    p.ctrl.EPS_disturb = real("control", "EPS_disturb");
    p.ctrl.V_Microblow0 = real("control", "V_Microblow0");
    if (table.contains("output")) {
        section("output", {"format", "physical_interval", "statistics_interval"});
        if (table["output"].as_table()->contains("physical_interval")) {
            p.physical_interval = integer("output", "physical_interval");
            if (p.physical_interval < 0)
                throw std::runtime_error("output.physical_interval must be non-negative");
        }
        if (table["output"].as_table()->contains("statistics_interval")) {
            p.statistics_interval = integer("output", "statistics_interval");
            if (p.statistics_interval <= 0)
                throw std::runtime_error("output.statistics_interval must be positive");
        }
        auto format = table["output"]["format"].value<std::string>();
        if (!format) throw std::runtime_error("output.format requires a string");
        if (lowercase(*format) == "ascii") p.output_format = OutputFormat::Ascii;
        else if (lowercase(*format) == "binary") p.output_format = OutputFormat::Binary;
        else if (lowercase(*format) == "both") p.output_format = OutputFormat::Both;
        else throw std::runtime_error("output.format must be ascii, binary, or both");
    }
    return p;
}

} // namespace

Params Params::load(const std::string& dir) {
    Params p{};
    const auto config = join_path(dir, "config.toml");
    if (std::filesystem::exists(config)) {
        p = load_toml(config);
    } else {

        {   // meshparm_channel_with_Microblow.txt —— 主程序 84-93 行
            const std::string f = join_path(dir, "meshparm_channel_with_Microblow.txt");
            LineCursor c(f);
            std::vector<std::string> t;
            c.skip(); t = c.tokens(); expect(f, t, 4);
            p.mesh.NX = to_int(f, t[0]); p.mesh.NY = to_int(f, t[1]);
            p.mesh.NZ = to_int(f, t[2]); p.mesh.INDEX_disturb_Z = to_int(f, t[3]);
            c.skip(); t = c.tokens(); expect(f, t, 3);
            p.mesh.NXS = to_int(f, t[0]); p.mesh.NYS = to_int(f, t[1]);
            p.mesh.NZS = to_int(f, t[2]);
            c.skip(); t = c.tokens(); expect(f, t, 3);
            p.mesh.NX2 = to_int(f, t[0]); p.mesh.NY2 = to_int(f, t[1]);
            p.mesh.NZ2 = to_int(f, t[2]);
            c.skip(); t = c.tokens(); expect(f, t, 1);
            p.mesh.INDEX_Serial = to_int(f, t[0]);
        }

        {   // parm_channel_with_Microblow.txt —— 主程序 181-194 行
            const std::string f = join_path(dir, "parm_channel_with_Microblow.txt");
            LineCursor c(f);
            std::vector<std::string> t;
            c.skip(); t = c.tokens(); expect(f, t, 3);
            p.parm.XL = to_double(f, t[0]); p.parm.ZL = to_double(f, t[1]);
            p.parm.THETA = to_double(f, t[2]);
            c.skip(); t = c.tokens(); expect(f, t, 2);
            p.parm.INDEX_GRAD_PRE = to_int(f, t[0]); p.parm.Re = to_double(f, t[1]);
            c.skip(); t = c.tokens(); expect(f, t, 1);
            p.parm.INDEX_UNSTEADY_BC = to_int(f, t[0]);
            c.skip(); t = c.tokens(); expect(f, t, 3);
            p.parm.G_pre = to_double(f, t[0]); p.parm.a_pre0 = to_double(f, t[1]);
            p.parm.b_pre0 = to_double(f, t[2]);
            c.skip();
            for (int i = 0; i < 3; ++i) {       // DO I=1,3 READ a_V(I),b_V(I)
                t = c.tokens(); expect(f, t, 2);
                p.parm.a_V[i] = to_double(f, t[0]); p.parm.b_V[i] = to_double(f, t[1]);
            }
        }

        {   // control_channel_with_Microblow.txt —— 主程序 205-216 行
            const std::string f = join_path(dir, "control_channel_with_Microblow.txt");
            LineCursor c(f);
            std::vector<std::string> t;
            c.skip(); t = c.tokens(); expect(f, t, 2);
            p.ctrl.INDEX_CONT = to_int(f, t[0]); p.ctrl.EPS_init = to_double(f, t[1]);
            c.skip(); t = c.tokens(); expect(f, t, 4);
            p.ctrl.NSTEPS = to_int(f, t[0]); p.ctrl.NDUMFRQ = to_int(f, t[1]);
            p.ctrl.NPRNFRQ = to_int(f, t[2]); p.ctrl.DT = to_double(f, t[3]);
            c.skip(); t = c.tokens(); expect(f, t, 3);
            p.ctrl.NCOUNT_disturb = to_int(f, t[0]); p.ctrl.NDISTURBFRQ = to_int(f, t[1]);
            p.ctrl.EPS_disturb = to_double(f, t[2]);
            c.skip(); t = c.tokens(); expect(f, t, 1);
            p.ctrl.V_Microblow0 = to_double(f, t[0]);
            c.skip(); t = c.tokens(); expect(f, t, 1);
            p.ctrl.Nstat = to_int(f, t[0]);
        }
    }

    // 派生物理参数
    const double pi = std::acos(-1.0);
    p.THETA_RAD = p.parm.THETA * pi / 180.0;
    if (p.parm.INDEX_GRAD_PRE == 1) {
        p.GRAD_PRE_X = -std::cos(p.THETA_RAD);
        p.GRAD_PRE_Z = -std::sin(p.THETA_RAD);
        p.QXmean = p.QZmean = 0.0;
    } else {
        p.QXmean = std::cos(p.THETA_RAD);
        p.QZmean = std::sin(p.THETA_RAD);
        p.GRAD_PRE_X = p.GRAD_PRE_Z = 0.0;
    }

    return p;
}
