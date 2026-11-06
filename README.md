# channel-dns

C++17 / MPI 不可压缩槽道流直接数值模拟（DNS）求解器，采用 Fourier–Chebyshev 谱方法和 FFTW，支持 MPI、进程内多线程、场输出与续算。

## 编译

依赖 C++17 编译器（已验证 GCC 13）、CMake ≥ 3.16、MPI、FFTW3 和 pkg-config。Ubuntu 24.04：

```sh
sudo apt-get install build-essential cmake pkg-config libopenmpi-dev openmpi-bin libfftw3-dev
./build.sh -j 8
```

生成 Release（`-O3`）程序 `build-channel-dns/bin/channel_dns`。默认使用线程池；若需每个 rank 在主线程直接计算，使用 `CHANNEL_DNS_PURE_MPI=ON ./build.sh -j 8`。

## 运行

在仓库根目录执行单进程示例（32³ 网格，5 步）：

```sh
mkdir -p runs/serial
cp examples/config.toml runs/serial/config.toml
./build-channel-dns/bin/channel_dns --timing runs/serial
```

MPI 双进程：

```sh
mkdir -p runs/mpi
cp examples/config.toml runs/mpi/config.toml
sed -i 's/INDEX_Serial = 1/INDEX_Serial = 0/' runs/mpi/config.toml
mpiexec -n 2 ./build-channel-dns/bin/channel_dns --threads 1 --timing runs/mpi
```

通过 `--threads N`（或 `-t N`）指定每 rank 的计算 worker 数，例如 `--threads 8`。未指定时使用 CPU 总逻辑线程数（检测失败回退为 1），另有主控制线程；MPI 多 rank 建议显式指定线程数并分配足够 CPU 核。纯 MPI 构建始终在主线程计算。运行目录保存 `field_S_*.dat` 续算场、`field_P_*.dat` 物理场、`EKT*.dat` 和 FFTW wisdom。修改 `config.toml` 的网格、时间步、输出频率等参数；设置 `INDEX_CONT=1` 从已有场续算，`output.format` 可选 `ascii`、`binary`、`both`。

[Reτ=180 槽道算例](examples/Re_tau180/Re_tau180.toml)：计算域 `4πh × 2h × 2πh`，基础网格 `128 × 129 × 128`，非线性乘积网格 `196 × 193 × 196`，固定压强梯度。

## 实测性能

纯 MPI 和线程池的布局测量、精度检查及复现命令见[性能记录](doc/performance.md)。

## 测试与文档

```sh
# 测试另需 Python ≥ 3.11、Git；首次配置会下载 googletest
./test.sh -j 8
```

- [数学算法：方程、空间离散、压力校正与时间推进](doc/algorithm.md)
- [C++ 设计：模块、计算流程、数据布局与访存优化](doc/design.md)
- [详细使用：配置、输出、续算、线程与 MPI 部署](doc/usage.md)
- [性能记录与测量口径](doc/performance.md)
