#pragma once

#include "fft/transforms.hpp"
#include "solver/fields.hpp"

namespace vector_ops {

// Reused by the timestepper; transient callers can omit the workspace argument.
struct ProductWorkspace {
    std::array<std::vector<double>, 3> a, b;
    std::vector<double> product;
    std::array<std::vector<double>, 2> cross_extra;
    void prepare(const Mesh &mesh);
};

void derivative(int axis, ScalarField &out, const ScalarField &in, const Mesh &mesh,
                const Transforms &transforms);
void gradient(VectorField &out, const ScalarField &in, const Mesh &mesh,
              const Transforms &transforms);
void divergence(ScalarField &out, const VectorField &in, const Mesh &mesh,
                const Transforms &transforms);
void curl(VectorField &out, const VectorField &in, const Mesh &mesh, const Transforms &transforms);
void enqueue_curl(VectorField &out, const VectorField &in, const Mesh &mesh,
                  const Transforms &transforms);
void laplacian(ScalarField &out, const ScalarField &in, const Mesh &mesh,
               const Transforms &transforms);
void divergence_boundary(ScalarBoundary &out, const VectorField &in, const Mesh &mesh,
                         const Transforms &transforms);
double max_distance(const ScalarField &lhs, const ScalarField &rhs, const Mesh &mesh);
double flow_rate(const ScalarField &in, const Mesh &mesh);
void dot_product(ScalarField &out, const VectorField &a, const VectorField &b, const Mesh &mesh,
                 Transforms &transforms, ProductWorkspace *workspace = nullptr);
void cross_product(VectorField &out, const VectorField &a, const VectorField &b, const Mesh &mesh,
                   Transforms &transforms, ProductWorkspace *workspace = nullptr);
void enqueue_cross_product(VectorField &out, const VectorField &a, const VectorField &b,
                           const Mesh &mesh, Transforms &transforms, ProductWorkspace &workspace);

} // namespace vector_ops
