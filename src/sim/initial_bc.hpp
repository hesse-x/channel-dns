#pragma once

#include <random>

#include "solver/fields.hpp"

class Mesh;
class Transforms;
struct Params;

namespace initial_bc {

using DisturbanceEngine = std::mt19937_64;

double centered_uniform(DisturbanceEngine& engine);

void initialize_velocity(VectorField& velocity, const Params& params,
                         const Mesh& mesh, Transforms& transforms);
void initialize_velocity(VectorField& velocity, const Params& params,
                         const Mesh& mesh, Transforms& transforms,
                         DisturbanceEngine& engine);
void disturb_velocity(VectorField& velocity, double amplitude, const Params& params,
                      const Mesh& mesh, Transforms& transforms, DisturbanceEngine& engine);
void velocity_boundary(VectorBoundary& boundary, double time,
                       const Params& params, const Mesh& mesh);

} // namespace initial_bc
