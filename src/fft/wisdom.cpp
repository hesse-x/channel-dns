#include "fft/wisdom.hpp"

#include <fftw3.h>

#include <filesystem>
#include <iostream>
#include <system_error>

namespace wisdom {

bool import(const std::string& file) {
    int ok = fftw_import_wisdom_from_filename(file.c_str());
    if (ok)
        std::cout << "channel-dns: FFTW wisdom loaded from " << file << "\n";
    else
        std::cout << "channel-dns: no FFTW wisdom loaded from " << file << "\n";
    return ok != 0;
}

void export_wisdom(const std::string& file) {
    const std::filesystem::path target(file);
    const std::filesystem::path temporary(file + ".tmp");
    const int ok = fftw_export_wisdom_to_filename(temporary.string().c_str());
    std::error_code error;
    if (ok) std::filesystem::rename(temporary, target, error);
    if (ok && !error) {
        std::cout << "channel-dns: FFTW wisdom saved to " << file << "\n";
    } else {
        std::filesystem::remove(temporary, error);
        std::cerr << "channel-dns: failed to export FFTW wisdom to " << file << "\n";
    }
}

}  // namespace wisdom
