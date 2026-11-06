// Mesh：网格尺寸、MPI 分块表与本 rank 坐标切片。
#pragma once

#include "core/params.hpp"

#include <iosfwd>
#include <vector>

class Mesh {
public:
    // XL/ZL 在 Parm 中（坐标计算需要），与 MeshParm 一并传入。
    Mesh(const MeshParm& mp, double XL, double ZL, int my_id, int nproc);

    // 网格尺寸
    int NX, NXS, NX2, NXH, NXS2P1, NXSP, NXSP2P1;
    int NY, NYS, NY2, NYP, NY2P;
    int NZ, NZS, NZ2, NZH, NZS2P1;
    int NXTNZ;
    int INDEX_Serial, my_id, NPROC, INDEX_MPI0, NPART;

    // 分块表（从零开始编号）
    std::vector<int> LSUBX, NXSPT, LSUBX2, NXSP2P1T;
    std::vector<int> LSUBY, NYPT, LSUBY2, NY2PT;

    // 坐标（本 rank 持有的切片）
    std::vector<double> XC, YC, ZC;   // XC(NX), YC(0:NYP) 本块, ZC(NZ)

    int x_begin(int rank) const;
    int x_count(int rank) const;
    int y_begin(int rank) const;
    int y_count(int rank) const;
    int y2_begin(int rank) const;
    int y2_count(int rank) const;
    bool owns_zero_x_mode() const;
    int zero_mode_owner() const;

    void print_tables(std::ostream& os) const;   // 验收用：打印全部分块表与尺寸
};
