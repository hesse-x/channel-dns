#pragma once

#include <iomanip>
#include <locale>
#include <sstream>
#include <string>

// Preserve the space-padded rank field (minimum width 2).
inline std::string rank_filename(const char* prefix, int rank, const char* suffix) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << prefix << std::setw(2) << rank << suffix;
    return out.str();
}
