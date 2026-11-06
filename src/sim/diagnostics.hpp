#pragma once

#include <filesystem>
#include "sim/statistics.hpp"
#include "solver/fields.hpp"

class Transforms;
class MpiContext;

struct DiagnosticValues {
    double ekt = 0.0;         // Unweighted sum(v*v) over grid nodes.
    double tke = 0.0;         // Volume average of 0.5*sum((u-<u>_xz)^2).
    double divergence_l2 = 0.0; // sqrt(volume average of (div u)^2).
    bool finite = true;
};

class Diagnostics {
public:
    Diagnostics(const Mesh& mesh, Transforms& transforms, const MpiContext* mpi);
    DiagnosticValues evaluate(const std::vector<PlaneMoments>& planes,
                              const VectorField& velocity, bool spectral_finite);
    void append(const std::filesystem::path& directory, std::int64_t step, double time,
                const DiagnosticValues& values) const;
private:
    const Mesh& m_;
    Transforms& t_;
    const MpiContext* mpi_;
    std::vector<double> weights_;
    ScalarField divergence_;
    PhysicalScalarField physical_divergence_;
};
