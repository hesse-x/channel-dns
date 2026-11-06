#pragma once

#include <cstdint>
#include <array>
#include <filesystem>
#include <vector>

#include "sim/physical_fields.hpp"

class Mesh;

// Plane sums: u,v,w,p; uu,vv,ww,uv,vw,wu; pu,pv,pw,pp; u^3,v^3,w^3,u^4,v^4,w^4.
struct PlaneMoments {
    std::array<double, 20> sum{};
    bool finite = true;
};

class Statistics {
public:
    explicit Statistics(const Mesh& mesh);

    void accumulate(const PhysicalVectorField& velocity,
                    const PhysicalScalarField& pressure, double dpdx, double dpdz);
    const std::vector<PlaneMoments>& reduce(const PhysicalVectorField& velocity,
                    const PhysicalScalarField& pressure, bool higher_moments, bool validate);
    void accumulate(const std::vector<PlaneMoments>& planes, double dpdx, double dpdz);
    void write(const std::filesystem::path& run_dir) const;

    std::int64_t count() const { return count_; }
    double dpdx_mean() const { return dpdx_mean_; }
    double dpdz_mean() const { return dpdz_mean_; }

private:
    const Mesh& m_;
    std::int64_t count_ = 0;
    double dpdx_mean_ = 0.0, dpdz_mean_ = 0.0;
    std::vector<double> velocity_mean_, pressure_mean_;
    std::vector<double> velocity_second_, pressure_velocity_, pressure_second_;
    std::vector<double> velocity_third_, velocity_fourth_;
    std::vector<PlaneMoments> plane_scratch_;
};
