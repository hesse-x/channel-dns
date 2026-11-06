# 使用说明

程序是 Fortran / 纯 MPI 求解器，默认使用单线程 FFTW，也可用 CMake 参数 `-DCHANNEL_DNS_FFT_BACKEND=FFTPACK` 选择 FFTPACK 内核；后者无需 FFTW 依赖。配置解析由 `src/toml.f90` 实现，无额外库依赖。

## 构建

```sh
./build.sh -j 8
# 输出：build-channel-dns/bin/channel_dns
./test.sh -B build-test -j 8
```

`FC` / `MPIFC` 可用于指定编译器，`MPIEXEC` 指定启动器。CMake 参数 `BUILD_LEGACY_REFERENCE=OFF` 可仅构建主求解器与相应测试。默认 legacy 基准需要 Python 3.11；运行测试和比较工具也需要 Python 3.11。

## 配置与运行

将 `examples/config.toml` 或 `examples/Re_tau180/Re_tau180.toml` 复制为运行目录中的 `config.toml`。配置表包括：

| 表 | 内容 |
| --- | --- |
| `mesh` | 粗/细网格 `NX/NY/NZ`、`NX2/NY2/NZ2`，谱截断 `NXS/NYS/NZS`，`INDEX_Serial` 和扰动开关 |
| `parm` | 周期长度、Re、流向角、固定压力梯度/固定流量选择、压力和速度边界系数 |
| `control` | 初场/续算、步数、时间步、输出频率、扰动及统计起始步数 |
| `output` | 输出格式、物理场与统计写出频率 |

`config.toml` 使用内置的配置子集解析器，语法范围如下：

- 单层表和裸键，例如 `[mesh]`、`NX = 32`；表名及键名大小写不敏感。
- 十进制整数及带 `e` / `E` 指数的实数，例如 `32`、`-0.5`、`3.0e-3`。
- 单行数值数组，例如 `a_V = [1, 1.0, 1e0]`；速度边界数组必须有三个元素。
- 单引号或双引号字符串，例如 `format = "binary"`，字符串中不使用转义。
- 空行、空格、制表符和 `#` 行尾注释；引号内的 `#` 属于字符串内容。

每行最多 4096 个字符。多行数组、嵌套表、数组表、带引号的键、字符串转义、布尔值、日期、非十进制整数和数字分隔符不在支持范围内，会明确报错。未知表、未知字段、重复名称、类型错误、整数越界和非有限参数均拒绝。不存在 `config.toml` 时读取三个文本参数文件，文件名见 README。

```sh
# 单 rank：INDEX_Serial = 1
./build-channel-dns/bin/channel_dns --timing runs/serial
# 多 rank：INDEX_Serial = 0
mpiexec --bind-to core --map-by core -n 4 ./build-channel-dns/bin/channel_dns --timing runs/mpi
```

rank 数不能超过 `min(NXS+1, NY+1, NY2+1)`。`--threads N` / `-t N` 接受线程数参数，但不产生计算线程。`--timing` 报告整个推进循环的最大 rank 墙钟耗时；`--benchmark-warmup N` 额外报告跳过 N 步后、不含场输出/诊断的计算阶段耗时，供比较脚本使用。

## 输出与续算

- `field_S_*.dat` / `.bin`：谱速度续算场。
- `field_P_*.dat` / `.bin`：物理速度和旋转形式压力。
- `diagnostics.csv`、`EKT*.dat`：能量、平面脉动能和散度诊断。
- `Result/0`～`Result/4`：压力梯度及平面一至四阶矩。
- `channel_dns.wisdom*`：各 rank FFTW 计划缓存。

rank 文件名保留宽度至少为 2 的空格填充，例如 `field_S_ 0.bin`。`output.format` 为 `ascii`、`binary` 或 `both`。二进制为小端 128 字节头加 binary64 数据。写出先完成临时文件，再替换目标文件。

`NDUMFRQ` 控制续算场写出，`NPRNFRQ` 控制诊断；`output.physical_interval` 默认跟随 `NPRNFRQ`，设为 0 可禁用物理场输出。`output.statistics_interval` 控制统计写出频率。统计在全局步数大于 `Nstat` 时每步采样。

设置 `INDEX_CONT = 1` 从同一运行目录续算，需要保持相同网格和 rank 数。文本模式优先 `.dat`；二进制和双格式模式优先 `.bin`。不要保留不同时间的混合格式检查点。二进制读取会校验网格、rank、载荷长度、数值及跨 rank 的步数和时间。

续算只保存速度，不保存多步历史或统计累积，因此每次续算重新启动时间推进系数和统计累积。随机扰动使用 Fortran 随机数生成器，若要对拍精度，必须使用共同的确定性续算场。

## legacy 精度与性能比较

按 README 中的 `tools/compare_legacy.py` 命令运行。工具统一读取/生成 TOML，再生成 legacy 参数文件；参考程序启动时以各自运行目录为工作目录。所有初始检查点使用相同精度的文本表示，工具记录 SHA-256，并分别保留日志和 JSON 结果。

参考程序使用 MPI 和 FFTPACK。`channel_dns_legacy_aligned` 仅用于比较相同零波数约束下的精度，性能基准使用 `channel_dns_legacy`。性能记录使用同一编译器、Release 优化、相同 MPI rank 数和绑核，交替重复运行。
