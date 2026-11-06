# channel-dns

C++17 / MPI 不可压缩槽道流直接数值模拟（DNS）求解器，采用 Fourier–Chebyshev 谱方法和 FFTW，支持 MPI、进程内多线程、场输出与续算。

## 纯 MPI 与混合模式

两种模式使用相同的数值算法，均支持单 rank / 多 rank 运行、场输出、统计与续算；区别在于 rank 内的计算执行方式。混合模式为 **MPI + C++ 线程池**。

| 项目 | 纯 MPI | 混合模式（默认） |
| --- | --- | --- |
| 编译选项 `CHANNEL_DNS_MPI_ONLY` | `ON` | `OFF` |
| rank 内计算 | 主线程直接执行，无计算 worker 和任务队列 | 每 rank 一个线程池，另有主线程负责控制和 MPI |
| 并行方式 | 通过增加 MPI rank 数并行 | MPI rank 间并行 + rank 内多线程，也可单 rank 多线程 |
| `--threads N` | 忽略，每 rank 主线程计算 | 指定每 rank 的计算 worker 数；`--threads 1` 仍使用线程池 |
| FFTW | 单线程 FFTW，不依赖 `fftw3_threads` | 需要 `fftw3_threads`，FFT 和主要计算复用线程池 |

这是编译期选择，不能通过 `--threads` 切换后端。两种模式的性能取决于网格、rank / 线程布局和 CPU 绑定，实测对比见[性能记录](doc/performance.md)。

部署时，**集群环境建议先尝试纯 MPI 模式**，通过多个 MPI rank 分配跨节点计算；**单机建议使用混合模式构建下的纯线程模式**，即单个 rank 配合 `--threads N`，并在 `config.toml` 中设置 `INDEX_Serial = 1`。例如 `./build-channel-dns/bin/channel_dns --threads 8 --timing /path/to/run` 使用 8 个计算 worker。最终可根据实际网格和硬件的测量结果调整布局。

## 编译

依赖 C++17 编译器（已验证 GCC 13）、CMake ≥ 3.16、MPI、FFTW3 和 pkg-config。Ubuntu 24.04：

```sh
sudo apt-get install build-essential cmake pkg-config libopenmpi-dev openmpi-bin libfftw3-dev

# 混合模式（默认）
./build.sh -j 8

# 纯 MPI
CHANNEL_DNS_MPI_ONLY=ON ./build.sh -j 8
```

上述命令均生成 Release（`-O3`）程序 `build-channel-dns/bin/channel_dns`，后一次构建会替换该目录中的程序；再次执行不带环境变量的 `./build.sh` 会恢复混合模式。`-j 8` 指编译并行任务数，运行时线程数由 `--threads` 指定。

如需同时保留两个版本，可使用独立 CMake 构建目录：

```sh
# 混合模式：build-hybrid/src/channel_dns
cmake -S . -B build-hybrid -DCMAKE_BUILD_TYPE=Release \
  -DCHANNEL_DNS_MPI_ONLY=OFF -DENABLE_UNIT_TESTS=OFF
cmake --build build-hybrid -j 8

# 纯 MPI：build-pure-mpi/src/channel_dns
cmake -S . -B build-pure-mpi -DCMAKE_BUILD_TYPE=Release \
  -DCHANNEL_DNS_MPI_ONLY=ON -DENABLE_UNIT_TESTS=OFF
cmake --build build-pure-mpi -j 8
```

Ubuntu 的 `libfftw3-dev` 同时提供单线程 FFTW 与 `fftw3_threads`；纯 MPI 构建只需前者。

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
