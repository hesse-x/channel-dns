#pragma once

#include <cstdint>
#include <filesystem>

#include "sim/physical_fields.hpp"
#include "core/params.hpp"
#include "solver/fields.hpp"

class Mesh;
class Transforms;
class MpiContext;

class FieldIO {
public:
    FieldIO(const Mesh& mesh, Transforms& transforms, std::filesystem::path run_dir,
            const MpiContext* mpi = nullptr, OutputFormat format = OutputFormat::Ascii);

    const PhysicalSnapshot& to_physical(const VectorField& velocity,
                                        const ScalarField& pressure, bool* finite = nullptr);
    void append_ekt(std::int64_t ntpit, double time, const PhysicalSnapshot& snapshot);
    void write_physical(std::int64_t ntpit, double time, const PhysicalSnapshot& snapshot);
    void write_restart(std::int64_t ntpit, double time, const VectorField& velocity);
    void read_restart(std::int64_t& ntpit, double& time, VectorField& velocity) const;

private:
    OutputFormat format_;
    const Mesh& m_;
    const MpiContext* mpi_;
    Transforms& transforms_;
    std::filesystem::path run_dir_;
    PhysicalSnapshot snapshot_;
    ScalarField spectral_workspace_;
};
