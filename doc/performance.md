# 性能测量记录

测试机器为 Intel Core i7-14700，28 个逻辑 CPU、20 个物理核，单 NUMA 节点；GCC 13.3.0、OpenMPI 4.1.6、FFTW 3.3.10。网格设置为 `NXS=NZS=grid/2`、`NYS=grid-2`，非线性网格参数等于 grid，`DT=0.003`。

## 纯 MPI 8 rank 与单 rank 8 线程

上述机器、Release `-O3`，对比两个独立构建：纯 MPI 使用 `CHANNEL_DNS_MPI_ONLY=ON`，8 rank、每 rank 主线程计算；线程池使用 `OFF`，1 rank、`--threads 8`，另有主控制线程。OpenMPI 均使用 `--bind-to core --map-by slot:PE=N`，纯 MPI 的 N=1，线程池 N=8，分配 8 个物理核；没有逐 worker 绑核。i7-14700 是异构 CPU，结果依赖绑核及机器负载。

使用相同零随机扰动配置、100 步、3 次测量取中位数，计时轮关闭场输出与统计。先零步预热生成 wisdom，每轮交替布局顺序；各网格依次运行，没有并发运行其他基准。以下是 `--timing` 的时间步循环耗时，取各 rank 最大值，包含首两步求解器 setup，排除 MPI 启动、构造、FFTW planning 与初始场准备；不包含完整进程启动和退出耗时。

| 网格参数 | 纯 MPI：8 rank | 线程池：1 rank × 8 worker | 纯 MPI / 线程池耗时 |
| --- | ---: | ---: | ---: |
| 32³ | 0.068824 s | 0.081173 s | 0.848× |
| 64³ | 0.502298 s | 0.558020 s | 0.900× |
| 128³ | 6.375658 s | 5.386789 s | 1.184× |

比值大于 1 表示线程池更快，小于 1 表示纯 MPI 更快。这里比较的是两种完整执行布局，包括它们不同的变换与通信路径，不能只归因于线程调度。

各轮原始循环耗时，单位秒：

| 网格参数 | 布局 | 第 1 次 | 第 2 次 | 第 3 次 |
| --- | --- | ---: | ---: | ---: |
| 32³ | 纯 MPI 8×1 | 0.073275 | 0.068824 | 0.064512 |
| 32³ | 线程池 1×8 | 0.079127 | 0.081173 | 0.086545 |
| 64³ | 纯 MPI 8×1 | 0.498957 | 0.511193 | 0.502298 |
| 64³ | 线程池 1×8 | 0.583521 | 0.558020 | 0.554572 |
| 128³ | 纯 MPI 8×1 | 6.138285 | 6.687946 | 6.375658 |
| 128³ | 线程池 1×8 | 5.346738 | 5.995272 | 5.386789 |

独立精度轮均推进 10 步后输出，与纯 MPI 8 rank 的场比较：三种网格的 u/v/w/P 最大绝对差均为 0；EKT 最大绝对差分别为 `1.125e-14`、`3.062e-13`、`5.728e-11`，以 `1e-10` 相对容差检查全部通过。并行归约顺序不同会导致 EKT 舍入差异。

复现（在仓库根目录执行，需要 gcc/MPI/FFTW；各工作目录须不存在）：

```sh
cmake -S . -B /tmp/channel-threaded-release -DCMAKE_BUILD_TYPE=Release \
  -DCHANNEL_DNS_MPI_ONLY=OFF -DENABLE_UNIT_TESTS=OFF
cmake --build /tmp/channel-threaded-release -j 8
cmake -S . -B /tmp/channel-pure-release -DCMAKE_BUILD_TYPE=Release \
  -DCHANNEL_DNS_MPI_ONLY=ON -DENABLE_UNIT_TESTS=OFF
cmake --build /tmp/channel-pure-release -j 8
for grid in 32 64 128; do
  python3 tests/tools/bench_hybrid.py /tmp/channel-threaded-release/src/channel_dns \
    --pure-mpi-executable /tmp/channel-pure-release/src/channel_dns \
    --grid "$grid" --steps 100 --repeats 3 --layouts 8x1,1x8 \
    --work-dir "/tmp/channel-layout-$grid"
done
```

`--pure-mpi-executable` 使每 rank 单线程布局使用纯 MPI 程序；不传此参数时，各布局均使用第一个参数指定的程序（默认线程池构建的 `8x1` 仍有每 rank 一个 worker）。脚本保存逐轮日志、CSV 和精度轮输出。它也支持 `--layouts 8x1,4x2,2x4,1x8` 及额外的 `--serial-baseline` 单线程基线。

这些零随机扰动算例没有验证一般湍流负载或跨节点扩展；非线性网格与粗网格尺寸相同，未额外配置完整 3/2 去混叠。见 [算法说明](algorithm.md)。
