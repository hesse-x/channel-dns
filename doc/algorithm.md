# 数学算法与离散求解

本文按当前 Fortran 实现说明算法。代码入口为 [`main.f90`](../src/main.f90)，线性求解在 [`dns_equations.f90`](../src/dns_equations.f90)，变换在 [`fft.f90`](../src/fft.f90)。模块、数据布局与 MPI 实现见 [设计说明](design.md)。以下标准不可压缩方程对应默认 `G_pre=0` 的配置；非零 `G_pre` 会改变压力算子。

## 1. 控制方程和边界

计算域为两个周期方向 x、z 和壁法向 $y\in[-1,1]$，周期长度分别为 `XL`、`ZL`。无量纲不可压缩 Navier–Stokes 方程写成旋转形式：

$$
\partial_t\boldsymbol u
=\boldsymbol u\times\boldsymbol\omega-\nabla P
+\frac{1}{Re}\Delta\boldsymbol u-\boldsymbol G,
\qquad
\boldsymbol\omega=\nabla\times\boldsymbol u,
\qquad \nabla\cdot\boldsymbol u=0.
$$

其中 $\boldsymbol u=(u,v,w)$，$\boldsymbol G=(G_x,0,G_z)$ 表示平均压力梯度，独立通过零模响应处理。利用

$$
(\boldsymbol u\cdot\nabla)\boldsymbol u
=\nabla(|\boldsymbol u|^2/2)-\boldsymbol u\times\boldsymbol\omega,
$$

可知代码压力变量 P 对应旋转形式中的修正压力：按该方程解释，$P=p+|\boldsymbol u|^2/2$。物理场输出直接变换求解器的 P，没有另行扣除动能项。

每个速度分量使用 `a_V[c]`、`b_V[c]` 给定的两壁 Robin 条件：

$$
a_c u_c(\pm1)+b_c\partial_n u_c(\pm1)=h_c^\pm,
\qquad \partial_n|_{+1}=\partial_y,\quad
\partial_n|_{-1}=-\partial_y.
$$

默认 $a_c=1,b_c=0$，即指定壁面速度。当前 `velocity_boundary()` 将两壁法向速度都设为 `V_Microblow0`，其余分量设为零；时间参数目前未用于构造时变边界。数组 wall=0 对应 $y=+1$，wall=1 对应 $y=-1$，不能按变量名把它们理解为几何下壁和上壁。

## 2. Fourier–Chebyshev 空间离散

场表示为截断展开：

$$
f(x,y,z)\approx
\sum_{k_x=-K_x}^{K_x}\sum_{k_z=-K_z}^{K_z}
\sum_{j=0}^{M}\widehat f_{k_x,j,k_z}
T_j(y)e^{i(\alpha k_x x+\beta k_z z)},
\quad \alpha=2\pi/XL,\ \beta=2\pi/ZL,
$$

其中 $K_x=\mathrm{NXS}$、$K_z=\mathrm{NZS}$、$M=\mathrm{NYS}$。利用实场的共轭对称性，只存非负 x 波数，z 存储 $-K_z,\ldots,K_z$。粗网格 x/z 等距，y 使用 Chebyshev–Gauss–Lobatto 点：

$$
y_\ell=\cos(\pi\ell/\mathrm{NY}),\quad \ell=0,\ldots,\mathrm{NY}.
$$

实际代码 x/z 坐标从一个网格步长开始，周期末端与原点等价。`NX/NY/NZ` 指物理变换网格参数；y 实际有 `NY+1` 个点；`NXS/NYS/NZS` 指保留模态，两者不应混淆。

### 导数与变换

Fourier 导数直接乘波数：$\partial_x\widehat f=i\alpha k_x\widehat f$，$\partial_z\widehat f=i\beta k_z\widehat f$。Chebyshev 一阶导数系数 d 由高阶向低阶递推：

$$
d_M=0,\quad d_{M-1}=2M\widehat f_M,
\quad d_j=d_{j+2}+2(j+1)\widehat f_{j+1},\quad
j=M-2,\ldots,0,
$$

最后将 $d_0$ 除以 2。递推对应 `dns_vector.f90` 及融合求解中的 y 导数。

FFTW 使用 x 向 r2c/c2r、z 向复数 FFT、y 向 `REDFT00`（DCT-I）。FFTPACK 后端对应使用 RFFT、CFFT 和 COST。两个后端的底层变换均不归一化，FFT 模块在展开/收谱时处理缩放及端点权重，最终谱系数符合上述直接求和约定。

### 非线性伪谱计算

每步先在谱空间计算旋度，再将速度和旋度的三个分量变换到 `NX2 × (NY2+1) × NZ2` 网格，逐点计算 $\boldsymbol N=\boldsymbol u\times\boldsymbol\omega$，随后正变换并截取保留模态。padding、乘积网格和截谱共同控制混叠。

代码不会自动按 3/2 规则扩大乘积网格。`examples/config.toml` 的乘积网格与粗网格相同，性能对照配置使用 3/2 倍乘积网格。实际 DNS 应根据保留模态和所需混叠控制选择乘积网格。

## 3. 时间推进：一阶启动、二阶半隐式多步法

设当前已有 $\boldsymbol u^n$ 及非线性项 $\boldsymbol N^n$，步长 $\Delta t=\mathrm{DT}$。默认压力算子下，代码求解：

$$
\gamma\boldsymbol u^{n+1}-\frac1{Re}\Delta\boldsymbol u^{n+1}
+\nabla P^{n+1}=\boldsymbol R^n,
$$

平均压力梯度留给零模校正。第一步以一阶隐式黏性、显式非线性启动：

$$
\gamma=\frac1{\Delta t},\qquad
\boldsymbol R^n=\frac{\boldsymbol u^n}{\Delta t}+\boldsymbol N^n.
$$

第二步起采用 BDF2 时间导数和二阶显式非线性外推：

$$
\gamma=\frac{3}{2\Delta t},\qquad
\boldsymbol R^n=\frac{2\boldsymbol u^n}{\Delta t}
-\frac{\boldsymbol u^{n-1}}{2\Delta t}
+2\boldsymbol N^n-\boldsymbol N^{n-1}.
$$

代码记号为 `a0,a1,b0,b1`，分别乘当前速度、历史速度、当前非线性项和历史非线性项，且 $\gamma=a_0+a_1$。`velocity_g=Re*gamma`。第一、二步分别重建线性求解系数、影响函数及压力梯度响应，之后步长和参数不变时复用。

本实现使用固定步长，没有自适应 CFL 选步。续算只恢复速度和时间，不恢复多步历史，因此每次续算都重新一阶启动。

## 4. 压力与速度的模态方程

对一个 Fourier 模态定义

$$
\kappa^2=(\alpha k_x)^2+(\beta k_z)^2,
\qquad L_k=D_y^2-\kappa^2.
$$

由离散动量方程的散度得到压力方程，随后求解三分量速度 Helmholtz 方程：

$$
L_k\widehat P=\widehat{\nabla\cdot\boldsymbol R},\qquad
(L_k-Re\gamma)\widehat u_c
=Re\bigl(\widehat{\partial_c P}-\widehat R_c\bigr).
$$

代码允许将压力方程推广成 $(L_k-\mathrm{G\_pre})\widehat P=\widehat{\nabla\cdot\boldsymbol R}$；这不是默认 Poisson 方程，非零配置的物理含义需由使用者明确。

### tau 线性求解

通用列求解的形式为 $(D_y^2-\kappa^2-g)f=r$，附两壁条件。实现对系数方程使用两次积分关系，生成连接 $j-2,j,j+2$ 的递推，偶数/奇数系数分别处理，再用两条边界条件闭合低阶自由度；没有为每个模态构造稠密矩阵求逆。

`WL/WD/WR` 是积分关系的系数，对 $j\ge2$：

$$
W_L(j)=\frac1{4j(j-1)},\quad
W_D(j)=-\frac1{2(j^2-1)},\quad
W_R(j)=\frac1{4j(j+1)},
$$

其中 $W_L(2)$ 再乘 2 以匹配零阶 Chebyshev 归一化。`build_wn_ad()` 预计算消元系数，`solve_tau_tile()` 在每组 8 个相邻模态上完成递推和边界闭合。

## 5. 压力边界、影响函数及 sigma 校正

压力并非只解一次 Poisson 后便结束。有限阶 tau 离散下，压力边界和截断高阶方程会影响速度的离散散度，代码分以下阶段校正。

### 非线性压力边界

`TimeStepper::enqueue_pressure_boundary()` 在两壁计算

$$
B^n=\boldsymbol n\cdot\left[
\boldsymbol N^n-\frac1{Re}\nabla\times\boldsymbol\omega^n\right].
$$

具体 wall=0 为 $N_y+(\partial_x\omega_z-\partial_z\omega_x)/Re$ 在 $y=+1$ 的值，wall=1 为其在 $y=-1$ 的相反数。Chebyshev 端点值通过系数求和与交错求和获得。首步使用 $B^n$，后续使用 $2B^n-B^{n-1}$。默认压力条件为 $\partial_n P=B$；配置还允许 `a_pre0*P+b_pre0*∂nP=B`。当前边界构造没有附加任意时变壁面速度的时间导数项。

### 影响函数：消除两壁散度

`setup()` 暂时对压力使用 Dirichlet 条件，预计算两组齐次响应 $(P^{(0)},\boldsymbol u^{(0)})$、$(P^{(1)},\boldsymbol u^{(1)})$：压力在两壁分别取 `(1,0)`、`(0,1)`，速度为齐次壁面条件，体 RHS 为零。随后恢复配置的压力边界。

每步先求特解 $(P^*,\boldsymbol u^*)$。对非零 Fourier 模态，记

$$
d_\pm^*=\nabla\cdot\boldsymbol u^*|_{y=\pm1},\qquad
A_{\pm,j}=\nabla\cdot\boldsymbol u^{(j)}|_{y=\pm1}.
$$

求解两阶影响矩阵

$$
A\begin{bmatrix}c_0\\c_1\end{bmatrix}
=-\begin{bmatrix}d_+^*\\d_-^*\end{bmatrix},
$$

并叠加 $P=P^*+c_0P^{(0)}+c_1P^{(1)}$、$\boldsymbol u=\boldsymbol u^*+c_0\boldsymbol u^{(0)}+c_1\boldsymbol u^{(1)}$。矩阵的行列式及归一化系数在初始化时预存。速度齐次响应不会改变已给定的速度壁面条件；压力边界在此作为构造特解的辅助条件，最终压力还包含影响函数校正，不能把初始压力 Robin 数据当成独立于不可压缩约束的最终物理条件。

端点散度无需完整逆变换：利用 $T_j(1)=1$、$T_j(-1)=(-1)^j$、$T'_j(1)=j^2$、$T'_j(-1)=(-1)^{j-1}j^2$ 直接求和。

### sigma：修正 tau 最高阶残差

影响函数处理两壁约束后，`build_b_problem()` 还提供压力/速度 B 响应 $(P^B,\boldsymbol u^B)$。其压力 RHS 系数为 $r_0=M-1$，$r_j=2[M-((j+1)\bmod2)]$（$1\le j<M$），$r_M=0$；速度体 RHS 和两壁数据为零，经同一耦合求解得到响应。

对影响函数校正后的场，设 $w=-(\kappa^2+Re\gamma)$，代码计算

$$
\sigma_0=\frac{-w\widehat v_M-Re\widehat R_{v,M}}{S_0},\qquad
\sigma_1=\frac{-w\widehat v_{M-1}+Re(2M\widehat P_M-\widehat R_{v,M-1})}{S_1},
$$

其中

$$
S_0=Re+w\widehat v_M^B,\qquad
S_1=Re+w\widehat v_{M-1}^B-2ReM\widehat P_M^B.
$$

按 $q=(j+1)\bmod2$ 更新

$$
\widehat P_j\mathrel{+}=\sigma_q\widehat P_j^B,\quad
(\widehat u_j,\widehat w_j)\mathrel{+}=\sigma_q(\widehat u_j^B,\widehat w_j^B),\quad
\widehat v_j\mathrel{+}=\sigma_{1-q}\widehat v_j^B.
$$

这一步针对 Chebyshev tau 截断产生的最高阶残差，法向速度与压力/切向速度使用不同奇偶分支。实际热路径将这些校正融合在局部模态缓冲中，最后一次性写回；这与“解一个压力增量，再全场做一次梯度投影”的实现不同。

### 全局零模与压力规范

$(k_x,k_z)=(0,0)$ 下，纯 Neumann 压力有常数自由度。`solve_tau_tile()` 对零模修正 RHS 的相容性，并将压力零阶系数固定为零；这是固定 Chebyshev 常数系数的规范，不等同于对最终压力做全场体积均值扣除。影响函数的零模仅保留独立壁面约束，B 问题的一个零模分母也单独处理。

`correct_zero_mode()` 最后将 $\widehat v_0$ 设为两壁法向数据的平均值，将该模态其余 y 系数置零，落实均匀穿透流的无散约束。

## 6. 平均压力梯度和恒流量控制

周期压力无法表示 x/z 方向的线性平均压降，因而另行保存 $G_x,G_z$。求解器预先计算单位压力梯度驱动的零模速度响应 $u_G,w_G$。流量采用归一化截面平均，对零模 Chebyshev 系数为

$$
Q(f)=\frac12\int_{-1}^{1}f(y)\,dy
=-\sum_{j\ \mathrm{even}}\frac{\widehat f_j}{j^2-1}.
$$

固定压力梯度模式直接使用配置导出的 $G_x,G_z$；恒流量模式求

$$
G_x=\frac{Q_x^{target}-Q(u^*)}{Q(u_G)},\qquad
G_z=\frac{Q_z^{target}-Q(w^*)}{Q(w_G)},
$$

再叠加 $u=u^*+G_xu_G$、$w=w^*+G_zw_G$。响应的 RHS 为 `Re`，对应控制方程中的 $-\boldsymbol G$ 符号。只有零模所属 rank 做此计算，多 rank 时广播梯度供统计使用。

## 7. 数值验证边界

FFT、导数、tau 求解、边界及时间步由解析解与单元测试验证；不同 MPI 分块通过场输出比较检查一致性。当前性能记录使用零随机扰动算例，不能替代湍流算例的网格收敛、时间收敛、散度误差和去混叠验证。算法行为以代码及实际配置为准。
