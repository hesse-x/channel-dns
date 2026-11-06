#pragma once

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace test_support {
class TempDirectory {
public:
    TempDirectory() {
        auto pattern = (std::filesystem::temp_directory_path() / "channel-test-XXXXXX").string();
        std::vector<char> buffer(pattern.begin(), pattern.end());
        buffer.push_back('\0');
        const auto result = ::mkdtemp(buffer.data());
        if (!result) throw std::runtime_error("cannot create test directory");
        path = result;
    }
    ~TempDirectory() {
        std::error_code ignored;
        std::filesystem::remove_all(path, ignored);
    }
    TempDirectory(const TempDirectory&) = delete;
    TempDirectory& operator=(const TempDirectory&) = delete;
    std::filesystem::path path;
};
inline std::string read(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read test fixture: " + path.string());
    return {std::istreambuf_iterator<char>(in), {}};
}
inline void write(const std::filesystem::path& path, const std::string& data) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out << data;
    if (!out) throw std::runtime_error("cannot write test fixture: " + path.string());
}
inline std::string replace(std::string text, const std::string& before, const std::string& after) {
    const auto pos = text.find(before);
    if (pos == std::string::npos) throw std::runtime_error("missing replacement fixture: " + before);
    text.replace(pos, before.size(), after);
    return text;
}
} // namespace test_support
