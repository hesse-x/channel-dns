#pragma once

#include <array>
#include <string>

struct MeshParm {          // meshparm_channel_with_Microblow.txt
    int NX, NY, NZ, INDEX_disturb_Z;
    int NXS, NYS, NZS;
    int NX2, NY2, NZ2;
    int INDEX_Serial;
};

struct Parm {              // parm_channel_with_Microblow.txt
    double XL, ZL, THETA;
    int    INDEX_GRAD_PRE;
    double Re;
    int    INDEX_UNSTEADY_BC;
    double G_pre, a_pre0, b_pre0;
    std::array<double, 3> a_V, b_V;   // DO I=1,3 READ a_V(I),b_V(I)
};

struct Control {           // control_channel_with_Microblow.txt
    int    INDEX_CONT;
    double EPS_init;
    int    NSTEPS, NDUMFRQ, NPRNFRQ;
    double DT;
    int    NCOUNT_disturb, NDISTURBFRQ;
    double EPS_disturb;
    double V_Microblow0;
    int    Nstat;
};

enum class OutputFormat { Ascii, Binary, Both };

struct Params {
    OutputFormat output_format = OutputFormat::Ascii;
    // -1 preserves the NPRNFRQ cadence; physical_interval=0 disables dumps.
    int physical_interval = -1;
    int statistics_interval = -1;
    MeshParm mesh;
    Parm     parm;
    Control  ctrl;

    // 派生物理参数
    double THETA_RAD;                  // THETA * pi / 180
    double GRAD_PRE_X, GRAD_PRE_Z;     // INDEX_GRAD_PRE==1 分支
    double QXmean, QZmean;             // else 分支

    // 优先读取 dir/config.toml；不存在时兼容旧的三个参数文件。
    static Params load(const std::string& dir);
};
