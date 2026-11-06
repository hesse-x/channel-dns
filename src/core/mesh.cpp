#include "core/mesh.hpp"

#include <cmath>
#include <ostream>
#include <stdexcept>

Mesh::Mesh(const MeshParm& mp, double XL, double ZL, int my_id_, int nproc_)
    : NX(mp.NX), NXS(mp.NXS), NX2(mp.NX2),
      NY(mp.NY), NYS(mp.NYS), NY2(mp.NY2),
      NZ(mp.NZ), NZS(mp.NZS), NZ2(mp.NZ2),
      INDEX_Serial(mp.INDEX_Serial), my_id(my_id_), NPROC(nproc_) {
    if (NPROC <= 0 || my_id < 0 || my_id >= NPROC)
        throw std::invalid_argument("invalid MPI rank or size");
    if (INDEX_Serial != 0 && INDEX_Serial != 1)
        throw std::invalid_argument("INDEX_Serial must be 0 or 1");
    // ---- 主程序 94-99 行 ----
    NXH = NX / 2;  NZH = NZ / 2;
    NXS2P1 = 2 * NXS + 1;  NZS2P1 = 2 * NZS + 1;

    INDEX_MPI0 = (1 - INDEX_Serial) * my_id;
    NPART = INDEX_Serial + (1 - INDEX_Serial) * NPROC;
    if (NPART <= 0 || NPART > NXS + 1 || NPART > NY + 1 || NPART > NY2 + 1)
        throw std::invalid_argument("MPI size exceeds an x/y decomposition extent");

    // ---- X 分块表（主程序 101-117 行）----
    // NXSPT 先保存元素数量，再转换为本 rank 的最大下标。
    LSUBX.resize(NPART); NXSPT.resize(NPART);
    LSUBX2.resize(NPART); NXSP2P1T.resize(NPART);

    NXSP = (NXS + 1) / NPART;
    std::fill(NXSPT.begin(), NXSPT.end(), NXSP);
    int K = (NXS + 1) % NPART;
    if (K > 0)
        for (int i = NPART - K; i <= NPART - 1; ++i)
            NXSPT[i] = NXSP + 1;
    LSUBX[0] = 0;
    for (int i = 1; i <= NPART - 1; ++i)
        LSUBX[i] = LSUBX[i - 1] + NXSPT[i - 1];
    for (int i = 0; i < NPART; ++i) LSUBX2[i] = LSUBX[i] + LSUBX[i];
    for (int i = 0; i < NPART; ++i) { NXSPT[i] = NXSPT[i] - 1; NXSP2P1T[i] = NXSPT[i] + NXSPT[i] + 1; }
    NXSP = NXSPT[INDEX_MPI0];  NXSP2P1 = NXSP2P1T[INDEX_MPI0];

    // ---- Y 分块表（主程序 119-134 行）----
    LSUBY.resize(NPART); NYPT.resize(NPART);

    NYP = (NY + 1) / NPART;
    std::fill(NYPT.begin(), NYPT.end(), NYP);
    K = (NY + 1) % NPART;
    if (K > 0)
        for (int i = NPART - K; i <= NPART - 1; ++i)
            NYPT[i] = NYP + 1;
    LSUBY[0] = 0;
    for (int i = 1; i <= NPART - 1; ++i)
        LSUBY[i] = LSUBY[i - 1] + NYPT[i - 1];
    for (int i = 0; i < NPART; ++i) NYPT[i] = NYPT[i] - 1;
    NYP = NYPT[INDEX_MPI0];

    // ---- NY2 分块表（主程序 136-151 行）----
    LSUBY2.resize(NPART); NY2PT.resize(NPART);

    NY2P = (NY2 + 1) / NPART;
    std::fill(NY2PT.begin(), NY2PT.end(), NY2P);
    K = (NY2 + 1) % NPART;
    if (K > 0)
        for (int i = NPART - K; i <= NPART - 1; ++i)
            NY2PT[i] = NY2P + 1;
    LSUBY2[0] = 0;
    for (int i = 1; i <= NPART - 1; ++i)
        LSUBY2[i] = LSUBY2[i - 1] + NY2PT[i - 1];
    for (int i = 0; i < NPART; ++i) NY2PT[i] = NY2PT[i] - 1;
    NY2P = NY2PT[INDEX_MPI0];

    auto validate = [&](const std::vector<int>& begin, const std::vector<int>& upper,
                        int extent, const char* name) {
        int next = 0;
        for (int rank = 0; rank < NPART; ++rank) {
            const int count = upper[rank] + 1;
            if (count <= 0 || begin[rank] != next)
                throw std::invalid_argument(std::string("invalid ") + name + " decomposition");
            next += count;
        }
        if (next != extent)
            throw std::invalid_argument(std::string("incomplete ") + name + " decomposition");
    };
    validate(LSUBX, NXSPT, NXS + 1, "x");
    validate(LSUBY, NYPT, NY + 1, "y");
    validate(LSUBY2, NY2PT, NY2 + 1, "y2");

    NXTNZ = NX * NZ;

    // ---- 坐标（SETUP 16-36 行）----
    const double pi = std::acos(-1.0);

    // 周期网格坐标为 XC[i] = (i+1)*DX
    XC.resize(NX);
    const double dx = XL / static_cast<double>(NX);
    for (int i = 0; i < NX; ++i)
        XC[i] = static_cast<double>(i + 1) * dx;

    // YC0TONY(0:NY)：Chebyshev Gauss-Lobatto 点，对称填负，中点为 0
    std::vector<double> YC0TONY(NY + 1);
    for (int j = 0; j <= NY / 2 - 1; ++j) {
        YC0TONY[j] = std::cos(pi * static_cast<double>(j) / static_cast<double>(NY));
        YC0TONY[NY - j] = -YC0TONY[j];
    }
    YC0TONY[NY / 2] = 0.0;

    // 本 rank 切片：YC(0:NYP) = YC0TONY(LSUBY(INDEX_MPI0) : +NYP)，共 NYP+1 个点
    YC.assign(YC0TONY.begin() + LSUBY[INDEX_MPI0],
              YC0TONY.begin() + LSUBY[INDEX_MPI0] + NYP + 1);

    ZC.resize(NZ);
    const double dz = ZL / static_cast<double>(NZ);
    for (int k = 0; k < NZ; ++k)
        ZC[k] = static_cast<double>(k + 1) * dz;
}

namespace {
int checked_rank(int rank, int size) {
    if (rank < 0 || rank >= size) throw std::out_of_range("mesh decomposition rank");
    return rank;
}
}

int Mesh::x_begin(int rank) const { return LSUBX[checked_rank(rank, NPART)]; }
int Mesh::x_count(int rank) const { return NXSPT[checked_rank(rank, NPART)] + 1; }
int Mesh::y_begin(int rank) const { return LSUBY[checked_rank(rank, NPART)]; }
int Mesh::y_count(int rank) const { return NYPT[checked_rank(rank, NPART)] + 1; }
int Mesh::y2_begin(int rank) const { return LSUBY2[checked_rank(rank, NPART)]; }
int Mesh::y2_count(int rank) const { return NY2PT[checked_rank(rank, NPART)] + 1; }
bool Mesh::owns_zero_x_mode() const { return x_begin(INDEX_MPI0) == 0; }
int Mesh::zero_mode_owner() const {
    for (int rank = 0; rank < NPART; ++rank)
        if (x_begin(rank) == 0) return rank;
    throw std::logic_error("x decomposition has no zero-mode owner");
}

void Mesh::print_tables(std::ostream& os) const {
    auto put = [&](const char* name, const std::vector<int>& v) {
        os << name;
        for (int x : v) os << ' ' << x;
        os << '\n';
    };
    os << "NX " << NX << " NXS " << NXS << " NX2 " << NX2
       << " NXH " << NXH << " NXS2P1 " << NXS2P1
       << " NXSP " << NXSP << " NXSP2P1 " << NXSP2P1 << '\n';
    os << "NY " << NY << " NYS " << NYS << " NY2 " << NY2
       << " NYP " << NYP << " NY2P " << NY2P << '\n';
    os << "NZ " << NZ << " NZS " << NZS << " NZ2 " << NZ2
       << " NZH " << NZH << " NZS2P1 " << NZS2P1 << " NXTNZ " << NXTNZ << '\n';
    os << "INDEX_Serial " << INDEX_Serial << " my_id " << my_id
       << " NPROC " << NPROC << " INDEX_MPI0 " << INDEX_MPI0
       << " NPART " << NPART << '\n';
    put("LSUBX", LSUBX);   put("NXSPT", NXSPT);
    put("LSUBX2", LSUBX2); put("NXSP2P1T", NXSP2P1T);
    put("LSUBY", LSUBY);   put("NYPT", NYPT);
    put("LSUBY2", LSUBY2); put("NY2PT", NY2PT);
}
