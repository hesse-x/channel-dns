# channel-dns

Fortran / 纯 MPI 不可压缩槽道流 DNS 求解器，采用 Fourier–Chebyshev 谱方法，默认使用单线程 FFTW，也支持 FFTPACK 后端。求解器源码位于 `src/`，运行参数使用 TOML 配置，参考求解器位于 `legacy/`。

## 编译与运行

默认 FFTW 构建依赖 GNU Fortran、CMake ≥ 3.16、MPI Fortran（`mpi_f08`）、FFTW3 和 pkg-config。参考目标、测试和比较工具需要 Python ≥ 3.11。已验证 GNU 13 / Open MPI / FFTW 3.3.10。

```sh
sudo apt-get install gfortran cmake python3 pkg-config libopenmpi-dev openmpi-bin libfftw3-dev
./build.sh -j 8
```

`-j` 只控制编译并行度。每个 MPI rank 在主线程执行数值计算。

### 新建算例与输出目录

命令最后的路径是运行目录：程序从中读取 `config.toml`，并将检查点、场文件、诊断与统计结果写到该目录。目录可使用相对路径或绝对路径。

以下命令从仓库根目录执行，新建一个 8 rank 算例：

```sh
mkdir -p runs/case01
cp examples/config.toml runs/case01/config.toml
```

编辑 `runs/case01/config.toml` 中已有的对应字段，保留其余配置：

```toml
# [mesh] 中
INDEX_Serial = 0

# [control] 中
INDEX_CONT = 0
NSTEPS = 1000
DT = 3.0e-3
NDUMFRQ = 100
NPRNFRQ = 100
```

`INDEX_CONT = 0` 表示生成初场；这里推进 1000 步，每 100 步保存检查点并写出诊断。正常运行命令：

```sh
mpiexec --bind-to core --map-by core -n 8 \
  ./build-channel-dns/bin/channel_dns runs/case01
```

保存终端日志：

```sh
mpiexec --bind-to core --map-by core -n 8 \
  ./build-channel-dns/bin/channel_dns runs/case01 \
  > runs/case01/run.log 2>&1
```

单 rank 运行时，将该目录配置中的 `INDEX_Serial` 改为 `1`，直接执行：

```sh
./build-channel-dns/bin/channel_dns runs/case01
```

`--timing` 是可选参数，用于结束时报告整个推进循环的最大 rank 墙钟耗时，包含循环内的诊断与输出；正常计算无需指定。

### 续算

保留运行目录中的 `config.toml` 和全部 rank 的 `field_S_*.dat` / `.bin` 检查点，将 `[control]` 中对应字段修改为：

```toml
INDEX_CONT = 1
NSTEPS = 2000
```

随后使用相同 rank 数重新运行。例如，继续上述 8 rank 算例：

```sh
mpiexec --bind-to core --map-by core -n 8 \
  ./build-channel-dns/bin/channel_dns runs/case01
```

`NSTEPS = 2000` 表示从最近保存的检查点再推进 2000 步。网格和 rank 数必须与检查点一致；检查点按 `NDUMFRQ` 间隔覆盖保存，程序中断后恢复的是最近保存的状态。续算恢复速度和时间，多步法与统计累积重新启动。

### 独立构建目录

独立 CMake 构建生成 `build-release/src/channel_dns`：

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DENABLE_UNIT_TESTS=OFF
cmake --build build-release -j 8
```

### 切换 FFT 后端

后端在编译时通过 `CHANNEL_DNS_FFT_BACKEND` 选择，支持 `FFTW`（默认）和 `FFTPACK`，不需要修改 TOML 配置。

显式构建 FFTW 后端，需要 FFTW3 开发包和 pkg-config：

```sh
cmake -S . -B build-fftw -DCMAKE_BUILD_TYPE=Release \
  -DCHANNEL_DNS_FFT_BACKEND=FFTW
cmake --build build-fftw -j 8
# 可执行文件：build-fftw/src/channel_dns
```

构建 FFTPACK 后端，使用仓库内的 FFTPACK 内核，无需安装 FFTW：

```sh
cmake -S . -B build-fftpack -DCMAKE_BUILD_TYPE=Release \
  -DCHANNEL_DNS_FFT_BACKEND=FFTPACK
cmake --build build-fftpack -j 8
# 可执行文件：build-fftpack/src/channel_dns
```

也可以对已有构建目录重新运行 CMake，显式设置上述选项后重新编译。CMake 会缓存后端选择；切回 FFTW 时需传入 `-DCHANNEL_DNS_FFT_BACKEND=FFTW`。运行命令和配置相同，只需使用对应目录的可执行文件。两个后端共用 MPI 转置、归一化、Nyquist 处理和其余求解器优化。

## 配置与结果

将 [examples/config.toml](examples/config.toml) 复制到运行目录作为 `config.toml`，按算例修改以下参数。Reτ=180 算例见 [Re_tau180.toml](examples/Re_tau180/Re_tau180.toml)。

| 配置表 | 配置项 | 含义与设置 |
| --- | --- | --- |
| `mesh` | `NX`、`NY`、`NZ` | 粗网格参数，实际物理点数为 `NX × (NY+1) × NZ` |
| `mesh` | `NXS`、`NYS`、`NZS` | Fourier–Chebyshev 谱截断阶数 |
| `mesh` | `NX2`、`NY2`、`NZ2` | 非线性乘积网格，需根据谱截断选择去混叠分辨率 |
| `mesh` | `INDEX_Serial` | 单 rank 设 `1`，多个 MPI rank 设 `0` |
| `mesh` | `INDEX_disturb_Z` | 随机扰动中 z 向速度分量的幅度倍数 |
| `parm` | `XL`、`ZL`、`Re` | x/z 周期长度及雷诺数，均为正数 |
| `parm` | `THETA`、`INDEX_GRAD_PRE` | 流向角（度）；`1` 为固定单位压力梯度，`0` 为固定单位平均流量 |
| `parm` | `G_pre` | 压力方程附加系数，通常为 `0` |
| `parm` | `a_pre0`、`b_pre0` | 压力 Robin 边界系数，默认 `0`、`1` 为法向导数条件 |
| `parm` | `a_V`、`b_V` | u/v/w 三分量的 Robin 边界系数；默认 `[1, 1, 1]`、`[0, 0, 0]` 为指定壁面速度 |
| `parm` | `INDEX_UNSTEADY_BC` | `1` 每步刷新壁面边界，`0` 使用初始化边界；当前边界值为常量 |
| `control` | `INDEX_CONT` | `0` 生成初场，`1` 读取运行目录中的谱速度检查点续算 |
| `control` | `NSTEPS`、`DT` | 本次推进步数和固定时间步长 |
| `control` | `NDUMFRQ`、`NPRNFRQ` | 检查点保存间隔、EKT 与诊断文件写出间隔，单位均为步 |
| `control` | `EPS_init` | 初始随机扰动幅度，`0` 禁用 |
| `control` | `NCOUNT_disturb`、`NDISTURBFRQ`、`EPS_disturb` | 推进中注入随机扰动的次数、间隔步数和幅度；次数设 `0` 禁用 |
| `control` | `V_Microblow0` | 两壁施加的恒定 y 向速度 |
| `control` | `Nstat` | 全局步数大于此值时，每步采样统计 |
| `output` | `format` | `"ascii"`、`"binary"` 或 `"both"`，控制场文件格式 |
| `output` | `physical_interval` | 物理场输出间隔，默认跟随 `NPRNFRQ`，`0` 禁用 |
| `output` | `statistics_interval` | 统计写出间隔，默认跟随 `NPRNFRQ`，必须为正数 |

例如，使用二进制检查点、关闭物理场输出，并每 1000 步写出统计：

```toml
[output]
format = "binary"
physical_interval = 0
statistics_interval = 1000
```

`output` 表可整体省略，默认使用文本格式。多 rank 运行前设置 `INDEX_Serial = 0`，rank 数不能超过 `min(NXS+1, NY+1, NY2+1)`。续算需保持相同网格和 rank 数；检查点保存速度与时间，续算时重新启动多步法和统计累积。

结果位于运行目录：`field_S_*.dat` / `.bin` 为谱速度检查点，`field_P_*.dat` / `.bin` 为物理速度与压力，`diagnostics.csv` 和 `EKT*.dat` 为诊断，`Result/` 为平面统计。输出频率与续算细节见[使用说明](doc/usage.md)。

`control.NPRNFRQ` 控制 EKT 和 `diagnostics.csv` 的追加写出间隔，例如 `NPRNFRQ = 100` 为每 100 个全局步写一次。终端也会在保存检查点和本次计算最后一步打印 EKT，但这些额外打印不会改变 EKT 文件的写出间隔。物理场和统计文件的写出频率分别通过 `output.physical_interval`、`output.statistics_interval` 设置。


## 单元测试与 legacy 对照

```sh
./test.sh -j 8
```

包含 Fortran 编写的配置读取、FFT、Nyquist、谱导数、制造解 Helmholtz 求解、零波数压力/法向速度约束、向量恒等式、网格分块和场文件往返测试，覆盖 1/2/3/4 rank，并检查融合求解器与分步求解器的等价性、FFT 对齐回退；另外运行端到端续算、损坏输入拒绝与 legacy 精度测试，无需下载测试框架。

参考求解器的源码和参数文件位于 `legacy/`，文件校验值见 `legacy/SOURCE_MANIFEST.json`。设置 `BUILD_LEGACY_REFERENCE=OFF` 可关闭参考程序及其测试；默认构建提供：

- `channel_dns_legacy_raw`：参考驱动与数值内核，需显式构建此目标。
- `channel_dns_legacy`：参考数值内核与计时驱动。
- `channel_dns_legacy_aligned`：零波数约束与主求解器一致的参考程序，用于压力精度对照。

```sh
cmake --build build-release --target channel_dns_legacy channel_dns_legacy_aligned -j 8
python3 tools/compare_legacy.py \
  --modern build-release/src/channel_dns \
  --legacy build-release/legacy/channel_dns_legacy \
  --aligned-legacy build-release/legacy/channel_dns_legacy_aligned \
  --work-dir /tmp/channel-dns-comparison \
  --grids 32 64 128 --ranks 1 4 8 --steps 200 --warmup 100 --repeats 3
```

比较脚本从同一 TOML 生成 legacy 输入，两边读取逐字相同的确定性续算场；精度分别报告速度、参考压力、非零水平波数压力以及约束对齐后的压力。性能以固定 rank 数、CPU 绑核下的计算阶段耗时为主，剔除预热和文件输出，另存循环及进程耗时。完整结果见[精度与性能报告](doc/performance.md)。

## 代码格式

安装 fprettify（推荐 0.3.7）后，运行 `./format.sh`，按 `.fprettify.rc` 统一格式化 `src/` 和 `tests/fortran/` 中的 Fortran 文件。可用 `FPRETTIFY` 指定工具路径。

## 文档

- [数学算法](doc/algorithm.md)
- [Fortran 模块与数据布局](doc/design.md)
- [配置、输出、续算与部署](doc/usage.md)
- [精度与性能实测](doc/performance.md)
- [Reτ=180 算例](examples/Re_tau180/Re_tau180.toml)
