#include "io/field_io.hpp"
#include "core/rank_filename.hpp"
#include "core/cpu_device.hpp"

#include <array>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <system_error>

#include "core/mesh.hpp"
#include "core/mpi_context.hpp"
#include "fft/transforms.hpp"

namespace {

constexpr std::array<unsigned char, 8> kMagic{{'C', 'H', 'D', 'N', 'S', 0, 0, 0}};
constexpr std::uint32_t kVersion = 1;
constexpr std::uint32_t kRestart = 1;
constexpr std::uint32_t kPhysical = 2;
constexpr std::uint64_t kHeaderBytes = 128;

void put_u32(std::ostream &out, std::uint32_t value) {
    std::array<char, 4> bytes{};
    for (unsigned i = 0; i < bytes.size(); ++i)
        bytes[i] = static_cast<char>(value >> (8 * i));
    out.write(bytes.data(), bytes.size());
}

void put_u64(std::ostream &out, std::uint64_t value) {
    std::array<char, 8> bytes{};
    for (unsigned i = 0; i < bytes.size(); ++i)
        bytes[i] = static_cast<char>(value >> (8 * i));
    out.write(bytes.data(), bytes.size());
}

void put_i64(std::ostream &out, std::int64_t value) {
    put_u64(out, static_cast<std::uint64_t>(value));
}

void put_double(std::ostream &out, double value) {
    static_assert(sizeof(double) == sizeof(std::uint64_t), "binary64 required");
    std::uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    put_u64(out, bits);
}

// Batch payload writes to avoid one ostream call per scalar on large fields.
void write_doubles(std::ostream &out, const double *values, std::size_t count) {
    const std::uint32_t endian_probe = 1;
    if (*reinterpret_cast<const unsigned char*>(&endian_probe) == 1)
        out.write(reinterpret_cast<const char *>(values), count * sizeof(double));
    else
        for (std::size_t i = 0; i < count; ++i)
            put_double(out, values[i]);
}

class BinaryPayload {
  public:
    explicit BinaryPayload(std::ostream &out) : out_(out) {}
    void append(double value) {
        buffer_[count_++] = value;
        if (count_ == buffer_.size())
            flush();
    }
    void flush() {
        write_doubles(out_, buffer_.data(), count_);
        count_ = 0;
    }

  private:
    std::ostream &out_;
    std::array<double, 8192> buffer_{};
    std::size_t count_ = 0;
};

std::uint32_t get_u32(std::istream &in) {
    std::array<unsigned char, 4> b{};
    in.read(reinterpret_cast<char *>(b.data()), b.size());
    if (!in)
        throw std::runtime_error("truncated field header");
    std::uint32_t value = 0;
    for (unsigned i = 0; i < b.size(); ++i)
        value |= static_cast<std::uint32_t>(b[i]) << (8 * i);
    return value;
}

std::uint64_t get_u64(std::istream &in) {
    std::array<unsigned char, 8> b{};
    in.read(reinterpret_cast<char *>(b.data()), b.size());
    if (!in)
        throw std::runtime_error("truncated field data");
    std::uint64_t value = 0;
    for (unsigned i = 0; i < b.size(); ++i)
        value |= static_cast<std::uint64_t>(b[i]) << (8 * i);
    return value;
}

std::int64_t get_i64(std::istream &in) { return static_cast<std::int64_t>(get_u64(in)); }

double get_double(std::istream &in) {
    const std::uint64_t bits = get_u64(in);
    double value;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

struct Header {
    std::uint32_t kind;
    std::int64_t ntpit;
    double time;
    std::array<std::int64_t, 12> dimensions;
};

void write_header(std::ostream &out, std::uint32_t kind, std::int64_t ntpit, double time,
                  const Mesh &m, std::int64_t components) {
    out.write(reinterpret_cast<const char *>(kMagic.data()), kMagic.size());
    put_u32(out, kVersion);
    put_u32(out, kind);
    put_i64(out, ntpit);
    put_double(out, time);
    put_i64(out, m.NX);
    put_i64(out, m.NYP + 1);
    put_i64(out, m.NZ);
    put_i64(out, components);
    put_i64(out, kind == kRestart ? m.NZS2P1 : 0);
    put_i64(out, kind == kRestart ? m.NXSP + 1 : 0);
    put_i64(out, m.my_id);
    put_i64(out, m.NPROC);
    put_i64(out, kind == kRestart ? m.x_begin(m.INDEX_MPI0) : 0);
    put_i64(out, kind == kRestart ? m.NXSP + 1 : m.NX);
    put_i64(out, kind == kRestart ? 0 : m.y_begin(m.INDEX_MPI0));
    put_i64(out, kind == kRestart ? m.NYS + 1 : m.NYP + 1);
}

Header read_header(std::istream &in) {
    std::array<unsigned char, 8> magic{};
    in.read(reinterpret_cast<char *>(magic.data()), magic.size());
    if (!in || magic != kMagic)
        throw std::runtime_error("invalid field magic");
    if (get_u32(in) != kVersion)
        throw std::runtime_error("unsupported field version");
    Header h{};
    h.kind = get_u32(in);
    h.ntpit = get_i64(in);
    h.time = get_double(in);
    for (auto &value : h.dimensions)
        value = get_i64(in);
    return h;
}

template <typename Writer> void atomic_write(const std::filesystem::path &path, Writer writer) {
    const auto temp = std::filesystem::path(path.string() + ".tmp");
    try {
        std::ofstream out(temp, std::ios::binary | std::ios::trunc);
        if (!out)
            throw std::runtime_error("cannot open temporary output: " + temp.string());
        writer(out);
        out.flush();
        if (!out)
            throw std::runtime_error("failed writing temporary output: " + temp.string());
        out.close();
        std::filesystem::rename(temp, path);
    } catch (...) {
        std::error_code ignored;
        std::filesystem::remove(temp, ignored);
        throw;
    }
}

std::uint64_t checked_product(std::initializer_list<std::uint64_t> factors) {
    std::uint64_t result = 1;
    for (const auto value : factors) {
        if (value && result > std::numeric_limits<std::uint64_t>::max() / value)
            throw std::runtime_error("field dimensions overflow");
        result *= value;
    }
    return result;
}

} // namespace

FieldIO::FieldIO(const Mesh &mesh, Transforms &transforms, std::filesystem::path run_dir,
                 const MpiContext *mpi, OutputFormat format)
    : format_(format), m_(mesh), mpi_(mpi), transforms_(transforms), run_dir_(std::move(run_dir)),
      spectral_workspace_(make_scalar_field(mesh)) {
    const std::size_t count = static_cast<std::size_t>(m_.NZ) * (m_.NYP + 1) * m_.NX;
    for (auto &component : snapshot_.velocity)
        component.resize(count);
    snapshot_.pressure.resize(count);
}

const PhysicalSnapshot &FieldIO::to_physical(const VectorField &velocity,
                                             const ScalarField &pressure, bool *finite) {
    std::vector<int> slab_valid(finite ? m_.NZS2P1 : 0, 1);
    for (int component = 0; component < 4; ++component) {
        cpu_device().for_each(-m_.NZS, m_.NZS + 1, [&](int iz) {
            bool valid = true;
            for (int iy = 0; iy <= m_.NYS; ++iy)
                for (int ix = 0; ix <= m_.NXSP; ++ix) {
                    const auto value = component < 3 ? velocity.at_z(iz, iy, ix, component)
                                                     : pressure.at_z(iz, iy, ix, 0);
                    spectral_workspace_.at_z(iz, iy, ix, 0) = value;
                    if (finite) valid &= std::isfinite(value.real()) && std::isfinite(value.imag());
                }
            if (finite) slab_valid[iz + m_.NZS] &= valid;
        });
        auto& output = component < 3 ? snapshot_.velocity[component] : snapshot_.pressure;
        transforms_.channel_xyzfft('B', output.data(),
                                   spectral_workspace_.data());
    }
    if (finite) *finite = std::all_of(slab_valid.begin(), slab_valid.end(), [](int v) { return v; });
    return snapshot_;
}

void FieldIO::append_ekt(std::int64_t, double time, const PhysicalSnapshot &snapshot) {
    double local_ekt = 0.0;
    for (const double v : snapshot.velocity[1])
        local_ekt += v * v;
    const double ekt =
        mpi_ ? mpi_->reduce_sum(local_ekt, m_.zero_mode_owner(), "EKT reduce") : local_ekt;
    if (!m_.owns_zero_x_mode())
        return;
    std::ostringstream line;
    line << std::scientific << std::setprecision(17) << time << ' ' << ekt << '\n';
    std::ofstream out(run_dir_ / rank_filename("EKT", m_.my_id, ".dat"), std::ios::app);
    if (!out)
        throw std::runtime_error("cannot append EKT output");
    out << line.str();
    out.flush();
    if (!out)
        throw std::runtime_error("failed writing EKT output");
}

void FieldIO::write_physical(std::int64_t ntpit, double time, const PhysicalSnapshot &snapshot) {
    if (format_ != OutputFormat::Binary) {
        const auto path = run_dir_ / rank_filename("field_P_", m_.my_id, ".dat");
        atomic_write(path, [&](std::ostream &out) {
            out << std::scientific << std::setprecision(17) << ntpit << ' ' << time << '\n';
            // Physical velocity layout: component is the fastest index.
            for (std::size_t i = 0; i < snapshot.pressure.size(); ++i)
                for (int c = 0; c < 3; ++c)
                    out << snapshot.velocity[c][i] << '\n';
            for (double value : snapshot.pressure)
                out << value << '\n';
        });
    }
    if (format_ == OutputFormat::Ascii)
        return;
    const auto path = run_dir_ / rank_filename("field_P_", m_.my_id, ".bin");
    atomic_write(path, [&](std::ostream &out) {
        write_header(out, kPhysical, ntpit, time, m_, 4);
        for (const auto &component : snapshot.velocity)
            write_doubles(out, component.data(), component.size());
        write_doubles(out, snapshot.pressure.data(), snapshot.pressure.size());
    });
}

void FieldIO::write_restart(std::int64_t ntpit, double time, const VectorField &velocity) {
    std::string error;
    try {
        if (format_ != OutputFormat::Binary) {
            const auto path = run_dir_ / rank_filename("field_S_", m_.my_id, ".dat");
            atomic_write(path, [&](std::ostream &out) {
                out << std::scientific << std::setprecision(17) << ntpit << ' ' << time << '\n';
                for (int iz = 0; iz < m_.NZS2P1; ++iz)
                    for (int iy = 0; iy <= m_.NYS; ++iy)
                        for (int ix = 0; ix <= m_.NXSP; ++ix)
                            for (int c = 0; c < 3; ++c) {
                                const auto v = velocity(iz, iy, ix, c);
                                out << '(' << v.real() << ',' << v.imag() << ")\n";
                            }
            });
        }
        if (format_ != OutputFormat::Ascii) {
            const auto path = run_dir_ / rank_filename("field_S_", m_.my_id, ".bin");
            atomic_write(path, [&](std::ostream &out) {
                write_header(out, kRestart, ntpit, time, m_, 3);
                BinaryPayload payload(out);
                for (int iz = 0; iz < m_.NZS2P1; ++iz)
                    for (int iy = 0; iy <= m_.NYS; ++iy)
                        for (int ix = 0; ix <= m_.NXSP; ++ix)
                            for (int component = 0; component < 3; ++component) {
                                const auto value = velocity(iz, iy, ix, component);
                                payload.append(value.real());
                                payload.append(value.imag());
                            }
                payload.flush();
            });
        }
    } catch (const std::exception &e) {
        error = e.what();
    }
    if (mpi_) {
        const int success = mpi_->allreduce_min(error.empty() ? 1 : 0, "restart write status");
        mpi_->barrier("restart commit barrier");
        if (!success)
            throw std::runtime_error(error.empty() ? "restart write failed on another rank"
                                                   : error);
    } else if (!error.empty()) {
        throw std::runtime_error(error);
    }
}

void FieldIO::read_restart(std::int64_t &ntpit, double &time, VectorField &velocity) const {
    auto path = run_dir_ / rank_filename("field_S_", m_.my_id, ".dat");
    const auto binary_path = run_dir_ / rank_filename("field_S_", m_.my_id, ".bin");
    if (!std::filesystem::exists(path) ||
        (format_ != OutputFormat::Ascii && std::filesystem::exists(binary_path)))
        path = binary_path;
    std::ifstream in;
    Header h{};
    VectorField candidate = make_vector_field(m_);
    std::string error;
    try {
        in.open(path, std::ios::binary);
        if (!in)
            throw std::runtime_error("cannot open restart: " + path.string());
        std::array<unsigned char, 8> probe{};
        in.read(reinterpret_cast<char *>(probe.data()), probe.size());
        in.clear();
        in.seekg(0);
        if (probe != kMagic) {
            auto punctuation = [&](char expected) {
                char c;
                if (!(in >> c) || c != expected)
                    throw std::runtime_error("invalid text complex syntax");
            };
            auto number = [&]() {
                in >> std::ws;
                std::string token;
                while (in && in.peek() != std::char_traits<char>::eof()) {
                    const char c = static_cast<char>(in.peek());
                    if (std::isspace(static_cast<unsigned char>(c)) || c == ')' || c == ',')
                        break;
                    token += static_cast<char>(in.get());
                }
                for (char &ch : token)
                    if (ch == 'd' || ch == 'D')
                        ch = 'e';
                std::size_t used = 0;
                const double value = std::stod(token, &used);
                if (used != token.size() || !std::isfinite(value))
                    throw std::runtime_error("invalid text restart number");
                return value;
            };
            // Read the integer separately to retain the full checkpoint step range.
            if (!(in >> h.ntpit) || h.ntpit < 0)
                throw std::runtime_error("invalid restart step");
            in >> std::ws;
            if (in.peek() == ',')
                in.get();
            h.time = number();
            for (int iz = 0; iz < m_.NZS2P1; ++iz)
                for (int iy = 0; iy <= m_.NYS; ++iy)
                    for (int ix = 0; ix <= m_.NXSP; ++ix)
                        for (int c = 0; c < 3; ++c) {
                            in >> std::ws;
                            if (in.peek() == ',')
                                in.get();
                            punctuation('(');
                            const double real = number();
                            punctuation(',');
                            const double imag = number();
                            punctuation(')');
                            candidate(iz, iy, ix, c) = {real, imag};
                        }
            char extra;
            while (in.get(extra))
                if (!std::isspace(static_cast<unsigned char>(extra)))
                    throw std::runtime_error("extra text restart data");
        } else {
            h = read_header(in);
            const std::array<std::int64_t, 12> expected{
                {m_.NX, m_.NYP + 1, m_.NZ, 3, m_.NZS2P1, m_.NXSP + 1, m_.my_id, m_.NPROC,
                 m_.x_begin(m_.INDEX_MPI0), m_.NXSP + 1, 0, m_.NYS + 1}};
            if (h.kind != kRestart || h.ntpit < 0 || !std::isfinite(h.time) ||
                h.dimensions != expected)
                throw std::runtime_error("restart header does not match this rank and mesh");
            const auto values = checked_product({static_cast<std::uint64_t>(m_.NZS2P1),
                                                 static_cast<std::uint64_t>(m_.NYS + 1),
                                                 static_cast<std::uint64_t>(m_.NXSP + 1), 3, 2});
            const auto file_size = std::filesystem::file_size(path);
            if (values > (std::numeric_limits<std::uint64_t>::max() - kHeaderBytes) / 8 ||
                file_size != kHeaderBytes + values * 8)
                throw std::runtime_error("restart payload length mismatch");

            for (int iz = 0; iz < m_.NZS2P1; ++iz)
                for (int iy = 0; iy <= m_.NYS; ++iy)
                    for (int ix = 0; ix <= m_.NXSP; ++ix)
                        for (int component = 0; component < 3; ++component) {
                            const double real = get_double(in), imag = get_double(in);
                            if (!std::isfinite(real) || !std::isfinite(imag))
                                throw std::runtime_error("restart contains non-finite data");
                            candidate(iz, iy, ix, component) = {real, imag};
                        }
        }
    } catch (const std::exception &e) {
        error = e.what();
    }
    if (mpi_) {
        const int success = mpi_->allreduce_min(error.empty() ? 1 : 0, "restart read status");
        if (!success)
            throw std::runtime_error(error.empty() ? "restart validation failed on another rank"
                                                   : error);
        const auto steps = mpi_->allgather_i64(h.ntpit, "restart step consistency");
        std::uint64_t time_bits = 0;
        std::memcpy(&time_bits, &h.time, sizeof(time_bits));
        const auto times = mpi_->allgather_u64(time_bits, "restart time consistency");
        for (int rank = 0; rank < m_.NPROC; ++rank) {
            if (steps[rank] != h.ntpit || times[rank] != time_bits ||
                m_.x_begin(rank) + m_.x_count(rank) > m_.NXS + 1)
                throw std::runtime_error("restart ranks do not form one consistent checkpoint");
        }
    } else if (!error.empty()) {
        throw std::runtime_error(error);
    }
    velocity = std::move(candidate);
    ntpit = h.ntpit;
    time = h.time;
}
