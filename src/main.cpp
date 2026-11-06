// 运行入口：MPI 初始化、参数一致性检查、分块、模拟与 fail-fast 收尾。
#include <mpi.h>

#include <charconv>
#include <chrono>
#include <string_view>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "core/mesh.hpp"
#include "core/mpi_context.hpp"
#include "core/params.hpp"
#include "core/cpu_device.hpp"
#include "fft/transforms.hpp"
#include "solver/eqsolve_pv.hpp"
#include "solver/timestepper.hpp"
#include "sim/simulation.hpp"

namespace {

std::uint64_t parameter_fingerprint(const std::filesystem::path& dir) {
    constexpr std::uint64_t offset = 1469598103934665603ULL;
    constexpr std::uint64_t prime = 1099511628211ULL;
    std::uint64_t hash = offset;
    const bool toml = std::filesystem::exists(dir / "config.toml");
    for (const char* name : {"config.toml", "meshparm_channel_with_Microblow.txt",
                             "parm_channel_with_Microblow.txt",
                             "control_channel_with_Microblow.txt"}) {
        if ((std::string(name) == "config.toml") != toml) continue;
        std::ifstream in(dir / name, std::ios::binary);
        if (!in) throw std::runtime_error(std::string("cannot fingerprint parameter file: ") + name);
        char c = 0;
        while (in.get(c)) {
            hash ^= static_cast<unsigned char>(c);
            hash *= prime;
        }
        hash ^= 0xff;
        hash *= prime;
    }
    return hash;
}

std::string fmt(double v) {
    std::ostringstream os;
    os << std::scientific << std::setprecision(16) << v;
    return os.str();
}

void dump_params(std::ostream& out, const Params& p, int rank) {
    const MeshParm& m = p.mesh;
    const Parm& pa = p.parm;
    const Control& c = p.ctrl;
    out << "DUMP rank " << rank << " meshparm "
              << m.NX << ' ' << m.NY << ' ' << m.NZ << ' ' << m.INDEX_disturb_Z << ' '
              << m.NXS << ' ' << m.NYS << ' ' << m.NZS << ' '
              << m.NX2 << ' ' << m.NY2 << ' ' << m.NZ2 << ' ' << m.INDEX_Serial << '\n';
    out << "DUMP rank " << rank << " parm1 "
              << fmt(pa.XL) << ' ' << fmt(pa.ZL) << ' ' << fmt(pa.THETA) << ' '
              << pa.INDEX_GRAD_PRE << ' ' << fmt(pa.Re) << ' ' << pa.INDEX_UNSTEADY_BC << '\n';
    out << "DUMP rank " << rank << " parm2 "
              << fmt(pa.G_pre) << ' ' << fmt(pa.a_pre0) << ' ' << fmt(pa.b_pre0) << '\n';
    out << "DUMP rank " << rank << " parm3";
    for (int i = 0; i < 3; ++i)
        out << ' ' << fmt(pa.a_V[i]) << ' ' << fmt(pa.b_V[i]);
    out << '\n';
    out << "DUMP rank " << rank << " ctrl "
              << c.INDEX_CONT << ' ' << fmt(c.EPS_init) << ' '
              << c.NSTEPS << ' ' << c.NDUMFRQ << ' ' << c.NPRNFRQ << ' ' << fmt(c.DT) << ' '
              << c.NCOUNT_disturb << ' ' << c.NDISTURBFRQ << ' ' << fmt(c.EPS_disturb) << ' '
              << fmt(c.V_Microblow0) << ' ' << c.Nstat << '\n';
    out << "DUMP rank " << rank << " derived " << fmt(p.THETA_RAD);
    if (p.parm.INDEX_GRAD_PRE == 1)
        out << ' ' << fmt(p.GRAD_PRE_X) << ' ' << fmt(p.GRAD_PRE_Z);
    else
        out << ' ' << fmt(p.QXmean) << ' ' << fmt(p.QZmean);
    out << '\n';
}

void dump_mesh(std::ostream& out, const Mesh& g, int rank) {
    out << "DUMP rank " << rank << " sizes "
              << g.NXH << ' ' << g.NZH << ' ' << g.NXS2P1 << ' ' << g.NZS2P1 << ' '
              << g.INDEX_MPI0 << ' ' << g.NPART << ' '
              << g.NXSP << ' ' << g.NXSP2P1 << ' ' << g.NYP << ' ' << g.NY2P << ' '
              << g.NXTNZ << '\n';
    auto puti = [&](const char* tag, const std::vector<int>& v) {
        out << "DUMP rank " << rank << ' ' << tag;
        for (int x : v) out << ' ' << x;
        out << '\n';
    };
    puti("LSUBX", g.LSUBX);   puti("NXSPT", g.NXSPT);
    puti("LSUBX2", g.LSUBX2); puti("NXSP2P1T", g.NXSP2P1T);
    puti("LSUBY", g.LSUBY);   puti("NYPT", g.NYPT);
    puti("LSUBY2", g.LSUBY2); puti("NY2PT", g.NY2PT);
    auto putd = [&](const char* tag, const std::vector<double>& v) {
        out << "DUMP rank " << rank << ' ' << tag;
        for (double x : v) out << ' ' << fmt(x);
        out << '\n';
    };
    putd("XC", g.XC);
    putd("YC", g.YC);
    putd("ZC", g.ZC);
}

} // namespace

int main(int argc, char** argv) {
    int provided = MPI_THREAD_SINGLE;
    if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) != MPI_SUCCESS) return 1;
    if (provided < MPI_THREAD_FUNNELED) {
        std::cerr << "MPI_THREAD_FUNNELED is required for hybrid execution\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MpiContext mpi(MPI_COMM_WORLD);
    const int rank = mpi.rank(), nproc = mpi.size();

    bool print_fft = false, print_timing = false;
    std::string dir;
    int requested_threads = 0;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--print-fft")
            print_fft = true;
        else if (std::string(argv[i]) == "--timing")
            print_timing = true;
        else if (std::string(argv[i]) == "--threads" || std::string(argv[i]) == "-t") {
            int value = 0;
            bool valid = i + 1 < argc;
            if (valid) {
                const std::string_view text(argv[++i]);
                const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
                valid = error == std::errc{} && end == text.data() + text.size() && value > 0;
            }
            if (!valid) {
                if (rank == 0) std::cerr << "--threads requires a positive integer\n";
                MPI_Finalize();
                return 2;
            }
            requested_threads = value;
        }
        else if (!argv[i][0] || argv[i][0] != '-') {
            if (!dir.empty()) {
                if (rank == 0) std::cerr << "usage: channel_dns [--print-fft] [--timing] [--threads N] <run-dir>\n";
                MPI_Finalize();
                return 2;
            }
            dir = argv[i];
        } else {
            if (rank == 0) std::cerr << "unknown option: " << argv[i] << '\n';
            MPI_Finalize();
            return 2;
        }
    }
    if (dir.empty() || !std::filesystem::is_directory(dir)) {
        if (rank == 0) std::cerr << "usage: channel_dns [--print-fft] [--timing] [--threads N] <existing-run-dir>\n";
        MPI_Finalize();
        return 2;
    }

    int result = 0;
    try {
        const int threads = cpu_device(requested_threads).threads();
#ifdef CHANNEL_DNS_MPI_ONLY
        if (rank == 0) std::cout << "channel-dns: " << nproc
                                 << " MPI rank(s), inline computation on each rank's main thread\n";
#else
        if (rank == 0) std::cout << "channel-dns: " << nproc << " MPI rank(s), "
                                 << threads << " compute worker(s) per rank (control thread excluded)\n";
#endif
        Params params = Params::load(dir);
        const auto fingerprints = mpi.allgather_u64(parameter_fingerprint(dir),
                                                     "parameter fingerprint");
        for (const auto value : fingerprints)
            if (value != fingerprints.front())
                throw std::runtime_error("parameter files differ between MPI ranks");
        if (params.mesh.INDEX_Serial == 1 && nproc != 1)
            throw std::runtime_error("INDEX_Serial=1 requires exactly one MPI rank");
        if (params.mesh.INDEX_Serial == 0 && nproc < 2)
            throw std::runtime_error("INDEX_Serial=0 requires at least two MPI ranks");
        if (params.ctrl.NDUMFRQ <= 0 || params.ctrl.NPRNFRQ <= 0)
            throw std::runtime_error("NDUMFRQ and NPRNFRQ must be positive");
        if (params.ctrl.DT <= 0.0) throw std::runtime_error("DT must be positive");
        if (params.ctrl.NCOUNT_disturb < 0)
            throw std::runtime_error("NCOUNT_disturb must be non-negative");
        if (params.ctrl.NCOUNT_disturb > 0 && params.ctrl.NDISTURBFRQ <= 0)
            throw std::runtime_error("NDISTURBFRQ must be positive when disturbances are enabled");
        if (params.parm.INDEX_UNSTEADY_BC != 0 && params.parm.INDEX_UNSTEADY_BC != 1)
            throw std::runtime_error("INDEX_UNSTEADY_BC must be 0 or 1");
        Mesh mesh(params.mesh, params.parm.XL, params.parm.ZL, rank, nproc);

        std::string wisdom_name = "channel_dns.wisdom";
        if (threads > 1) wisdom_name += ".threads" + std::to_string(threads);
        if (nproc > 1)
            wisdom_name += ".np" + std::to_string(nproc) + ".rank" + std::to_string(rank);
        const auto wisdom = (std::filesystem::path(dir) / wisdom_name).string();
        Transforms fft(mesh, params.parm.XL, params.parm.ZL, params.parm.Re,
                       wisdom, &mpi);

        if (print_fft && rank == 0) {
            std::ostringstream buf;
            dump_params(buf, params, rank);
            dump_mesh(buf, mesh, rank);
            auto putd = [&](const char* tag, const std::vector<double>& v) {
                buf << "DUMP rank " << rank << ' ' << tag;
                for (double x : v) buf << ' ' << fmt(x);
                buf << '\n';
            };
            putd("WAVX", fft.WAVX);
            buf << "DUMP rank " << rank << " WAVZ";  // 存储偏移 zoff_，按 -NZS..NZS 打印
            for (int k = -mesh.NZS; k <= mesh.NZS; ++k)
                buf << ' ' << fmt(fft.WAVZ[k + fft.zoff()]);
            buf << '\n';
            putd("WL", fft.WL);
            putd("WD", fft.WD);
            putd("WR", fft.WR);
            std::cout << buf.str();
        }
        Simulation simulation(params, mesh, fft, dir, &mpi);
        simulation.initialize();
        if (print_timing) mpi.barrier("timestep timing start");
        const auto start = std::chrono::steady_clock::now();
        simulation.run();
        if (print_timing) {
            const double local = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
            double maximum = 0.0;
            if (MPI_Reduce(&local, &maximum, 1, MPI_DOUBLE, MPI_MAX, 0, mpi.comm()) != MPI_SUCCESS)
                throw std::runtime_error("timestep timing reduction failed");
            if (rank == 0) std::cout << "channel-dns: timestep wall time = "
                                     << std::setprecision(12) << maximum << " s\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "rank " << rank << ": " << e.what() << '\n';
        if (nproc > 1) {
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }
        result = 1;
    }

    MPI_Finalize();
    return result;
}
