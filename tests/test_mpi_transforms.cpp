#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/cpu_device.hpp"
#include "core/mesh.hpp"
#include "core/mpi_context.hpp"
#include "core/params.hpp"
#include "fft/transforms.hpp"

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

MeshParm mesh_parameters() {
    MeshParm p{};
    p.NX = 12;
    p.NY = 8;
    p.NZ = 12;
    p.NXS = 4;
    p.NYS = 6;
    p.NZS = 3;
    p.NX2 = 18;
    p.NY2 = 12;
    p.NZ2 = 18;
    p.INDEX_Serial = 0;
    return p;
}

double analytic(int ix, int global_y, int iz, int nx, int ny, int nz) {
    double x = 0.0;
    for (int mode = 0; mode <= 4; ++mode)
        x += static_cast<double>(mode + 1) *
             std::cos(2.0 * kPi * mode * static_cast<double>(ix) / nx);
    return x * std::cos(2.0 * kPi * static_cast<double>(global_y) / ny) *
           std::cos(2.0 * kPi * static_cast<double>(iz) / nz);
}

void check_transform(Transforms &transforms, const Mesh &mesh, bool fine) {
    const int nx = fine ? mesh.NX2 : mesh.NX;
    const int ny = fine ? mesh.NY2 : mesh.NY;
    const int nz = fine ? mesh.NZ2 : mesh.NZ;
    const int local_y = fine ? mesh.y2_count(mesh.my_id) : mesh.y_count(mesh.my_id);
    const int y_begin = fine ? mesh.y2_begin(mesh.my_id) : mesh.y_begin(mesh.my_id);
    std::vector<double> physical(static_cast<std::size_t>(nz) * local_y * nx);
    for (int iz = 0; iz < nz; ++iz)
        for (int iy = 0; iy < local_y; ++iy)
            for (int ix = 0; ix < nx; ++ix)
                physical[(static_cast<std::size_t>(iz) * local_y + iy) * nx + ix] =
                    analytic(ix, y_begin + iy, iz, nx, ny, nz);
    const auto original = physical;
    std::vector<std::complex<double>> spectral(static_cast<std::size_t>(mesh.NZS2P1) *
                                               (mesh.NYS + 1) * mesh.x_count(mesh.my_id));
    if (fine)
        transforms.channel_xyzfft2('F', physical.data(), spectral.data());
    else
        transforms.channel_xyzfft('F', physical.data(), spectral.data());

    double spectral_error = 0.0;
    for (int iz = -mesh.NZS; iz <= mesh.NZS; ++iz)
        for (int iy = 0; iy <= mesh.NYS; ++iy)
            for (int ix = 0; ix < mesh.x_count(mesh.my_id); ++ix) {
                const int global_x = mesh.x_begin(mesh.my_id) + ix;
                const double expected =
                    (iy == 2 && std::abs(iz) == 1)
                        ? (global_x == 0 ? 0.5 : 0.25 * static_cast<double>(global_x + 1))
                        : 0.0;
                const auto index = (static_cast<std::size_t>(iz + mesh.NZS) * (mesh.NYS + 1) + iy) *
                                       mesh.x_count(mesh.my_id) +
                                   ix;
                spectral_error = std::max(spectral_error, std::abs(spectral[index] - expected));
            }
    if (spectral_error > 2e-13)
        throw std::runtime_error(
            std::string(fine ? "fine" : "coarse") +
            " distributed forward transform mismatch: " + std::to_string(spectral_error));

    if (fine)
        transforms.channel_xyzfft2('B', physical.data(), spectral.data());
    else
        transforms.channel_xyzfft('B', physical.data(), spectral.data());
    double roundtrip_error = 0.0;
    for (std::size_t i = 0; i < physical.size(); ++i)
        roundtrip_error = std::max(roundtrip_error, std::abs(physical[i] - original[i]));
    if (roundtrip_error > 2e-12)
        throw std::runtime_error(
            std::string(fine ? "fine" : "coarse") +
            " distributed roundtrip mismatch: " + std::to_string(roundtrip_error));
    const std::complex<double> sentinel(17.0, -23.0);
    std::vector<std::complex<double>> interleaved(spectral.size() * 3, sentinel);
    if (fine)
        transforms.channel_xyzfft2('F', physical.data(), interleaved.data() + 1, 3);
    else
        transforms.channel_xyzfft('F', physical.data(), interleaved.data() + 1, 3);
    for (std::size_t i = 0; i < spectral.size(); ++i) {
        if (std::abs(interleaved[3 * i + 1] - spectral[i]) > 2e-13 ||
            interleaved[3 * i] != sentinel || interleaved[3 * i + 2] != sentinel)
            throw std::runtime_error("strided forward transform corrupted vector components");
    }
    const auto saved = interleaved;
    if (fine)
        transforms.channel_xyzfft2('B', physical.data(), interleaved.data() + 1, 3);
    else
        transforms.channel_xyzfft('B', physical.data(), interleaved.data() + 1, 3);
    if (interleaved != saved)
        throw std::runtime_error("strided inverse mutated spectral input");
    for (std::size_t i = 0; i < physical.size(); ++i)
        if (std::abs(physical[i] - original[i]) > 2e-12)
            throw std::runtime_error("strided roundtrip mismatch");
}

void check_xz_transform(Transforms &transforms, const Mesh &mesh) {
    std::vector<double> physical(static_cast<std::size_t>(mesh.NZ) * mesh.NX);
    for (int iz = 0; iz < mesh.NZ; ++iz)
        for (int ix = 0; ix < mesh.NX; ++ix) {
            double x = 0.0;
            for (int mode = 0; mode <= mesh.NXS; ++mode)
                x += static_cast<double>(mode + 1) *
                     std::cos(2.0 * kPi * mode * static_cast<double>(ix) / mesh.NX);
            physical[static_cast<std::size_t>(iz) * mesh.NX + ix] =
                x * std::cos(2.0 * kPi * static_cast<double>(iz) / mesh.NZ);
        }
    std::vector<std::complex<double>> spectral(static_cast<std::size_t>(mesh.NZS2P1) *
                                               mesh.x_count(mesh.my_id));
    transforms.xzfft(physical.data(), spectral.data());
    double error = 0.0;
    for (int iz = -mesh.NZS; iz <= mesh.NZS; ++iz)
        for (int ix = 0; ix < mesh.x_count(mesh.my_id); ++ix) {
            const int global_x = mesh.x_begin(mesh.my_id) + ix;
            const double expected =
                std::abs(iz) == 1 ? (global_x == 0 ? 0.5 : 0.25 * static_cast<double>(global_x + 1))
                                  : 0.0;
            const auto index =
                static_cast<std::size_t>(iz + mesh.NZS) * mesh.x_count(mesh.my_id) + ix;
            error = std::max(error, std::abs(spectral[index] - expected));
        }
    if (error > 2e-13)
        throw std::runtime_error("distributed xz transform mismatch: " + std::to_string(error));
}

} // namespace

int main(int argc, char **argv) {
    if (argc > 1) cpu_device(std::stoi(argv[1]));
    else cpu_device(1);
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED)
        MPI_Abort(MPI_COMM_WORLD, 1);
    MpiContext mpi(MPI_COMM_WORLD);
    std::cout << "rank: " << mpi.rank() << "\tsize: " << mpi.size() << "\n";
    int result = 0;
    try {
        if (mpi.size() < 2 || mpi.size() > 4)
            throw std::runtime_error("MPI transform test requires 2, 3, or 4 ranks (MPI size=" +
                                     std::to_string(mpi.size()) + ", PMI_SIZE=" +
                                     (std::getenv("PMI_SIZE") ? std::getenv("PMI_SIZE") : "unset") +
                                     ")");
        Mesh mesh(mesh_parameters(), 2.0 * kPi, 2.0 * kPi, mpi.rank(), mpi.size());
        if (mesh.x_count(mpi.rank()) <= 0 || mesh.y_count(mpi.rank()) <= 0 ||
            mesh.y2_count(mpi.rank()) <= 0)
            throw std::runtime_error("empty decomposition block");
        if (mesh.owns_zero_x_mode() != (mpi.rank() == mesh.zero_mode_owner()))
            throw std::runtime_error("zero-mode ownership mismatch");
        if (!mesh.owns_zero_x_mode() && mesh.x_begin(mpi.rank()) == 0)
            throw std::runtime_error("non-owner starts at zero x mode");
        Transforms transforms(mesh, 2.0 * kPi, 2.0 * kPi, 100.0, "", &mpi);
        check_xz_transform(transforms, mesh);
        check_transform(transforms, mesh, false);
        check_transform(transforms, mesh, true);
        mpi.barrier("MPI transform test completion");
        if (mpi.rank() == 0)
            std::cout << "MPI transform test passed with " << mpi.size() << " ranks\n";
    } catch (const std::exception &e) {
        std::cerr << "rank " << mpi.rank() << ": " << e.what() << '\n';
        result = 1;
    }
    const int global_result = mpi.allreduce_min(result == 0 ? 1 : 0, "MPI test result");
    MPI_Finalize();
    return global_result ? 0 : 1;
}
