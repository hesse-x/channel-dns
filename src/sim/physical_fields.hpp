#pragma once

#include <array>
#include <cstddef>
#include <vector>

using PhysicalScalarField = std::vector<double>;
using PhysicalVectorField = std::array<PhysicalScalarField, 3>;

struct PhysicalSnapshot {
    PhysicalVectorField velocity;
    PhysicalScalarField pressure;
};

inline std::size_t physical_index(int iz, int iy, int ix, int ny_points, int nx) {
    return (static_cast<std::size_t>(iz) * ny_points + iy) * nx + ix;
}
