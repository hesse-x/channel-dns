# 详细使用说明

C++ / MPI 槽道流 DNS 求解器。使用 FFTW 进行变换，采用 TOML 配置，支持文本和二进制场文件。

## 1. 安装依赖

构建 C++ 程序需要：

- 支持 C++17 的编译器；当前已验证 GCC 13。
- CMake ≥ 3.16。
- MPI 开发库及运行工具。
- 系统 FFTW3 双精度开发库，以及 `pkg-config`。

本文命令均在仓库根目录执行。以下以 Ubuntu 24.04 为例：

```sh
sudo apt-get update
sudo apt-get install -y build-essential cmake pkg-config libopenmpi-dev openmpi-bin libfftw3-dev
```

运行全部测试还需要 Python ≥ 3.11 和 Git：

```sh
sudo apt-get install -y python3 git
```

FFTW 直接使用系统安装的库，无需 FFTW 子模块或 OCaml。TOML 解析器 toml++ 已随源码提供，无需额外安装。

## 2. 编译 Release 程序

在仓库根目录执行：

```sh
./build.sh
```

默认使用 4 个并行任务，以 **Release / `-O3`** 构建。也可以指定并行度：

```sh
./build.sh -j 8
```

脚本自动创建 `build-channel-dns/`，生成的程序为：

```text
build-channel-dns/bin/channel_dns
```

此脚本只构建 C++ 求解器，不构建测试。需要指定 MPI 实现时，可设置环境变量，例如：

```sh
MPICXX=/usr/bin/mpicxx.openmpi MPIEXEC=/usr/bin/mpiexec.openmpi ./build.sh
```

也可以直接使用 CMake，得到相同的程序路径：

```sh
cmake -S . -B build-channel-dns \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_UNIT_TESTS=OFF \
  -DCMAKE_RUNTIME_OUTPUT_DIRECTORY="$PWD/build-channel-dns/bin"
cmake --build build-channel-dns --target channel_dns --parallel 4
```

## 3. 跑一个单进程例子

以下命令均在仓库根目录执行。先为本次运行创建目录，复制完整的配置示例，再启动程序：

```sh
mkdir -p runs/serial
cp examples/config.toml runs/serial/config.toml
./build-channel-dns/bin/channel_dns runs/serial
```

[examples/config.toml](../examples/config.toml) 默认设置：

| 配置 | 值 | 含义 |
| --- | --- | --- |
| `mesh.INDEX_Serial` | `1` | 单进程运行 |
| `mesh.NX / NY / NZ` | `32 / 32 / 32` | 物理网格参数 |
| `control.INDEX_CONT` | `0` | 从初始场开始 |
| `control.NSTEPS` | `5` | 推进 5 步 |
| `control.DT` | `0.003` | 时间步长 |
| `control.NDUMFRQ` | `1` | 每步保存续算场 |
| `control.NPRNFRQ` | `1` | 每步记录 EKT、湍动能和散度残差 |
| `output.format` | `"ascii"` | 输出文本 |
| `output.physical_interval` | 未指定 | 默认跟随 `NPRNFRQ`；设为 `0` 禁止物理场输出，正整数指定独立间隔 |
| `output.statistics_interval` | 未指定 | 统计场覆盖保存间隔，默认跟随 `NPRNFRQ`；不影响逐步采样 |

正常运行时，终端按诊断间隔显示步数、时间、EKT、TKE、散度残差和有限值标志，最后时间约为 `0.015`。运行目录中会生成：

```text
runs/serial/
├── config.toml
├── channel_dns.wisdom
├── EKT 0.dat
├── diagnostics.csv
├── field_P_ 0.dat
└── field_S_ 0.dat
```

- `field_S_ 0.dat`：最新一次输出的谱空间速度场，可用于续算。
- `field_P_ 0.dat`：最新一次输出的物理空间速度和压力。
- `EKT 0.dat`：追加记录时间与壁法向速度平方和。
- `channel_dns.wisdom`：FFTW 规划缓存，后续运行会尝试复用。

文件名中的 rank 按最小两字符宽度补空格，操作文件时需要加引号，例如：

```sh
head -n 3 'runs/serial/field_S_ 0.dat'
cat 'runs/serial/EKT 0.dat'
```

场文件在每次输出时覆盖，EKT 追加写入。重复跑独立算例时，建议使用新的运行目录。

`diagnostics.csv` 同样追加保存，每行包括 `step,time,ekt_sum_v2,tke_plane_fluctuation_volume_mean,divergence_l2,finite`：

- EKT 沿用原程序定义，为全网格节点的 `sum(v²)`，没有体积归一化。
- TKE 为相对当前 x-z 平面平均速度的脉动动能体积平均 `0.5*<sum((u-<u>_xz)²)>_V`；它与相对长期平均场的湍动能定义不同。x-z 使用均匀平均，y 使用 Clenshaw–Curtis 求积。
- `divergence_l2` 为 `sqrt(<(div u)²>_V)`，衡量连续性误差，不是动量方程残差或稳态收敛指标。
- `finite` 同时检查输入谱系数、物理速度/压力、矩归约和散度，避免仅检查几个标量而漏掉被实变换丢弃的非法虚部。检查复用谱复制和矩计算；发生异常时先写诊断再终止。续算保存时及最后一步也执行检查。它不能保证发现两次检查之间出现又消失的异常。

Reτ=180 示例每 1000 步记录标量，每 10000 步覆盖保存续算场和统计场，不输出三维物理场：

```toml
[output]
format = "binary"
physical_interval = 0
statistics_interval = 10000
```

## 4. 跑一个 MPI 双进程例子

复制示例，将 `INDEX_Serial` 改为 `0`，使用两个 MPI 进程运行：

```sh
mkdir -p runs/mpi
cp examples/config.toml runs/mpi/config.toml
sed -i 's/INDEX_Serial = 1/INDEX_Serial = 0/' runs/mpi/config.toml
mpirun -np 2 ./build-channel-dns/bin/channel_dns --threads 1 runs/mpi
```

每个 rank 分别输出自己的场文件；EKT 由持有零模的 rank 输出，FFTW 缓存也按进程分别保存。`INDEX_Serial = 1` 要求恰好一个进程；`INDEX_Serial = 0` 要求至少两个进程，且分块不能为空。

可以对照以上两个例子的结果，比较工具会自动重组分布式物理场：

```sh
python3 tests/tools/compare_fields.py runs/serial runs/mpi
```

## 5. 修改配置、输出和续算

配置文件必须命名为运行目录下的 `config.toml`。`[mesh]`、`[parm]`、`[control]` 定义网格、物理参数和时间推进控制，所有这些参数都必须填写；`a_V` 和 `b_V` 是各含三个数字的数组。浮点指数使用 TOML 的 `e`，不能使用 `D`。

表名、参数名及 `output.format` 的选项值不区分 ASCII 大小写；仅大小写不同的重复定义会报错。程序拒绝未知键、缺失参数、错误类型及非有限数值。若 `config.toml` 不存在，仍支持读取原来的三个文件：

```text
meshparm_channel_with_Microblow.txt
parm_channel_with_Microblow.txt
control_channel_with_Microblow.txt
```

两种配置同时存在时，以 TOML 为准。MPI 各进程会检查配置文件内容是否一致。

### 输出格式

修改现有 `[output]` 表中的 `format`：

```toml
[output]
format = "both" # 可选 "ascii"、"binary"、"both"
```

| 模式 | 场文件 |
| --- | --- |
| `ascii`（默认） | `.dat` 文本，速度分量连续排列，复数表示为 `(实部,虚部)` |
| `binary` | `.bin` 二进制，包含版本、尺寸和分块信息，批量写出 |
| `both` | 同时输出 `.dat` 和 `.bin` |

EKT 和统计结果始终使用文本。文本输出采用足够精度保证 binary64 往返。

### 统计结果

累计步数 **大于** `control.Nstat` 时，每步采样；按 `output.statistics_interval` 覆盖保存，最后一步也保存。采样不依赖物理场或诊断输出频率。默认示例 `Nstat = 50000`，因此跑 5 步不会产生 `Result/`。若要在小例子中查看统计结果，将现有配置中的 `Nstat` 改为 `0`：

```sh
sed -i 's/Nstat = 50000/Nstat = 0/' runs/serial/config.toml
```

后续运行会在 `Result/0/` 保存样本数和平均压力梯度，在 `Result/1/` 至 `Result/4/` 保存均值及二至四阶矩。

### 从已有场继续运行

保留同一运行目录的续算文件，将 `INDEX_CONT` 改为 `1`，再启动程序：

```sh
sed -i 's/INDEX_CONT = 0/INDEX_CONT = 1/' runs/serial/config.toml
./build-channel-dns/bin/channel_dns runs/serial
```

`NSTEPS` 表示本次追加运行的步数；读取到的累计步数和时间会继续增加。网格与 MPI 进程数应保持一致。续算仅恢复速度场、步数和时间，并重新启动多步法历史和统计累计。

读取支持文本、新 `.bin` 文件，以及旧 C++ 二进制 `.dat` 文件。`ascii` 模式优先读取 `.dat`，不存在时读取 `.bin`；`binary` 和 `both` 模式优先读取 `.bin`，不存在时读取 `.dat`。切换输出模式不会删除另一种格式的旧文件，应避免混用不同时间的续算文件。文本没有网格元数据，只能校验数据量、数值及各 rank 的步数和时间；二进制还会校验网格和分块信息。

## 6. 运行全部测试

```sh
./test.sh
# 或指定并行度：
./test.sh -j 8
```

脚本在独立的 `build-test/` 中以 Debug 构建，启用断言，再运行全部 CTest 测例，包括：

- TOML 配置读取、类型、必填项、非法输入和优先级。
- 数组、网格、场布局、FFT、向量运算、求解器、时间推进及边界接口。
- 输出格式、损坏续算文件处理、统计数值和模拟生命周期。
- 1–4 进程 MPI 通信接口和分布式变换。
- 文本与二进制读写、续算一致性及小网格续算精度。

首次配置测试会从 GitHub 获取 googletest。离线环境可提供已有源码：

```sh
GTEST_SOURCE_DIR=/path/to/googletest ./test.sh
```

可用 `./test.sh -B /path/to/build` 指定测试构建目录；`MPICXX` 和 `MPIEXEC` 可指定 MPI 工具。失败日志位于测试构建目录下的 `Testing/Temporary/`。两个脚本均支持 `-h`，从其他当前目录调用时仍能定位仓库。

## 二进制格式 v1

文件头为 128 字节，对应 Python `struct` 格式 `<8sIIqd12q`。所有整数为小端，浮点使用小端 IEEE binary64。

| 偏移 | 类型 | 含义 |
| --- | --- | --- |
| 0 | 8 字节 | `CHDNS\0\0\0` |
| 8 | uint32 | 版本号 `1` |
| 12 | uint32 | 类型：`1` 续算场，`2` 物理场 |
| 16 | int64 | 累计步数 |
| 24 | binary64 | 时间 |
| 32 | 12 × int64 | NX、本地物理 y 点数、NZ、分量数、z 模态数、本地 x 模态数、rank、进程数、x 起点、x 数量、y 起点、y 数量 |

续算载荷的循环顺序为 z、y、本地 x、分量，分量变化最快；每个复数依次保存实部和虚部。物理场载荷依次保存三个连续的速度分量数组和压力数组，每个数组中 x 变化最快，其次本地 y、z。元数据的起点从零计数；续算场的 y 起点及数量描述完整谱空间 y 轴，物理场则描述本 rank 的物理空间分块。

## 纯 MPI 编译：主线程直接计算

默认构建使用线程池。开启 `CHANNEL_DNS_PURE_MPI` 后，每个 rank 的 `CpuDevice` 都在
主线程直接执行，头文件中的 `enqueue`、范围遍历和 kernel 可被编译器内联；没有任务队列、
worker 或提交时的 callable 类型擦除。`copy_n` / `fill_n` 直接调用标准算法，`sync()` 为空操作。
此构建不编译线程池实现，也不依赖 `fftw3_threads`，只使用单线程 FFTW。

```sh
CHANNEL_DNS_PURE_MPI=ON ./build.sh -j 8
```

也可以单独维护一个 CMake 构建目录：

```sh
cmake -S . -B build-pure-mpi -DCMAKE_BUILD_TYPE=Release \
  -DCHANNEL_DNS_PURE_MPI=ON -DENABLE_UNIT_TESTS=OFF
cmake --build build-pure-mpi -j 8
mpiexec -n 8 --bind-to core build-pure-mpi/src/channel_dns --timing /path/to/run
```

多 rank 配置使用 `INDEX_Serial = 0`；单 rank 也可使用此构建。
`--threads N` 在纯 MPI 构建中被忽略，`device.threads()` 固定返回 1，
表示当前 rank 的主线程。调用返回时结果已就绪；异常在调用处直接抛出。
这是编译期后端选择，不根据运行时线程数切换。重新用 `./build.sh` 构建会恢复默认线程池后端。
纯 MPI 回归可使用 `CHANNEL_DNS_PURE_MPI=ON ./test.sh -B build-test-pure-mpi`。

## MPI 与单例线程池混合并行

以下描述默认的线程池构建。每个 MPI rank 使用一个进程内单例线程池。`--threads N`（或 `-t N`）指定池内计算 worker 的数量；未指定时取 CPU 总逻辑线程数（`std::thread::hardware_concurrency()`，检测失败回退为 1），每 rank 独立使用该数量。设置为 8 时创建 8 个持久计算 worker，另有独立的主线程负责控制和 MPI。主线程不计入这个数量，也不领取线程池任务；`sync()` 只等待任务完成。MPI 只在主线程执行，程序请求 `MPI_THREAD_FUNNELED`。

执行层只公开 `CpuDevice`，线程池是其私有实现。计算代码通过 Device 提交、同步、遍历、复制和填充，不直接包含或访问线程池：

```cpp
auto& device = cpu_device();
device.enqueue(task_a);
device.enqueue(task_b);  // 等 task_a 完成后开始；主线程提交后立即返回
device.enqueue_for_each(0, count, [](int i) { /* 本批次内的独立工作并行执行 */ });
device.sync();          // 主线程等待；任务全部由内部 worker 执行
// 此后可在主线程读取结果；依赖前面批次的任务也可以提前保序入队
```

提交按 FIFO 保序：前一次 `enqueue` 或 `enqueue_for_each` 全部完成后，后一次提交才开始。一次范围提交内部的分块由多个 worker 并行执行，不保证索引完成顺序。提交与主线程异步，主线程可继续处理其他工作；`sync()` 阻塞主线程直到队列中的工作完成。异步任务引用的数据必须保留到 `sync()` 完成。任务异常在队列排空后由 `sync()` 传回主线程，失败不会取消后续提交。嵌套的同步循环在当前 worker 中执行，不创建第二个线程池；任务内部不提交新队列任务，也不调用 `sync()`。

范围提交以批次持有各分块的独立 callable，队列按批次分发分块，避免为每个分块单独包装 `std::function`。所有 callable 副本构造成功后才发布批次；复制失败时，该批次不会有部分任务进入队列。前一批任务持有的捕获对象释放后，才开始下一批。

`device.launch(Grid1D, kernel)` 提交分块 kernel；`device.for_each`、`device.copy_n`、`device.fill_n` 是同步便捷操作。Device 的公开头文件不包含线程池、互斥量、条件变量或 worker 定义，调度实现集中在 `src/core/cpu_device.cpp`。

单 rank 的细网格 FFT 将 x/z/y 各阶段作为保序批次提前入队，各 worker 使用单线程 FFTW plan 处理独立平面；plan 在初始化时创建并复用，执行时匹配输入/输出缓冲的对齐类别。多 rank 的 x 方向批量变换仍通过 FFTW 官方回调接口复用线程池，MPI 交换前由主线程等待本地计算完成。谱空间求解和主要向量运算也使用这一线程池。参见 [FFTW 自定义线程后端](https://www.fftw.org/fftw3_doc/Usage-of-Multi_002dthreaded-FFTW.html)。构建需要系统 `fftw3_threads` 库，Ubuntu 的 `libfftw3-dev` 包提供该库。

单 rank 时间步会提前提交 curl、非线性变换、叉积、RHS、求解、校正及缓冲交换，由队列保证依赖顺序；主线程在步末读取结果或输出前等待。单 rank 的流量校正不调用 MPI 广播。多 rank 路径仍保留通信边界的显式等待。

单机纯线程运行，配置使用 `INDEX_Serial = 1`：

```bash
./build-channel-dns/bin/channel_dns --threads 8 --timing /path/to/run
```

4 台机器、每台 8 个物理核的 OpenMPI 启动示例，配置使用 `INDEX_Serial = 0`：

```bash
mpiexec -n 4 --hostfile hosts \
  --map-by ppr:1:node:PE=8 --bind-to core --report-bindings \
  /shared/channel_dns --threads 8 --timing /shared/run
```

`hosts` 为四台节点的 hostfile，每节点至少提供 8 个核的 slots。上述映射为每节点分配一个 rank、每 rank 绑定 8 个核；工作线程继承该 rank 的 CPU 允许范围，由操作系统在范围内调度，没有额外的逐线程绑核。配置文件和可执行程序须在各节点相同路径可见，输出目录需可写。参见 [OpenMPI 启动与绑定说明](https://www.open-mpi.org/doc/v4.1/man1/mpirun.1.php)。

`--timing` 报告所有 rank 中最大的完整时间步循环耗时，包含循环内的计算、通信和已配置的输出，排除 FFTW planning、求解器构造及初始场准备。多线程 wisdom 按线程数另存。

按固定 8 核预算对比默认线程池构建的不同 rank/worker 布局（`8x1` 仍是每 rank 一个 worker）：

```bash
python3 tests/tools/bench_hybrid.py build-openmpi-release/src/channel_dns \
  --grid 128 --steps 50 --repeats 3 \
  --layouts 8x1,4x2,2x4,1x8 --work-dir /tmp/channel-hybrid-run
```

纯 MPI 构建只应使用每 rank 单线程的 `Nx1` 布局。工作目录必须不存在。脚本预热后依次测量各布局，保留每轮日志与 CSV，并比较不同布局的场和 EKT；计时轮不写场文件。该脚本用于单机同核数对照，不模拟跨节点网络，也不保证混合并行比纯 MPI 更快。

### CPU kernel 与 grid 接口

`CpuDevice` 是共享主机内存的 CPU 执行设备，内部引用同一个单例线程池，不创建第二组线程。`Grid1D` 描述逻辑工作项总数和每个块包含的工作项数，kernel 收到该块的起止索引。逻辑块可以多于工作线程；当前后端将相邻逻辑块合并为较大的队列任务，避免逐元素排队。

```cpp
#include "core/cpu_device.hpp"

auto& device = cpu_device();
device.launch(Grid1D{count, 4096}, [=](Block1D block) {
    for (auto i = block.begin; i < block.end; ++i) {
        const auto value = input[i];
        scratch[i] = value * scale;       // 数据准备
        output[i] = scratch[i] + bias;    // 同一块内继续计算并写回
    }
});
// 还可以提交其他独立 kernel。
device.sync();
```

不同 launch 按 FIFO 批次保序，前一批完成后才开始后一批；主线程读取结果前需要 `sync()`。块内是普通 CPU 顺序循环，没有 GPU 式线程组屏障或独立显存。当前提供一维 grid，更多维度可在线性工作项中映射。

梯度、散度和旋度已按谱空间 z 平面组织成 kernel，在块内直接读取原场、完成导数和组合并写回，减少分量提取和中间场。求解器按相邻 8 个模态使用 worker 局部缓冲，从速度 RHS 在局部计算散度，再完成压力、三分量速度与校正，只写回最终场；非线性物理缓冲跨时间步复用，FFT 直接读写交错向量的指定分量，省去谱分量中间场。历史场按生命周期交换缓冲；单 rank FFT 复用 y/zx 工作区，归一化与收谱融合，展开时只清零 padding。物理叉积一次读取六个输入，同时计算三个输出分量。大块填充、复制、物理场准备和统计计算也使用线程池。随机数序列、单列递推和共享文件流写入仍保持各自所需的顺序。

基准脚本加 `--serial-baseline` 可额外测量单 rank 单线程，用作加速比基线；其余布局仍须使用相同总核数。
