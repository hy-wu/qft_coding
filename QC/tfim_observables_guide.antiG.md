# 转动场伊辛模型（TFIM）的量子计算模拟：物理可观测量与算法原理

一维转动场伊辛模型（Transverse Field Ising Model, TFIM）是凝聚态物理和量子场论中最重要的基石模型之一。它不仅展示了丰富的量子相变（Quantum Phase Transition, QPT），还是量子计算动力学模拟（Quantum Dynamics Simulation）的经典基准。

本指南将深入拆解在量子计算机上模拟 TFIM 的**核心物理可观测量**及其物理意义，并详细阐述底层的**算法实现与模拟原理**。

---

## 一、 物理可观测量（Observables）与物理意义

在 TFIM 模拟中，体系的哈密顿量为：
$$H = -J \sum_{i=1}^{N-1} Z_i Z_{i+1} - g \sum_{i=1}^{N} X_i$$
其中 $Z_i, X_i$ 是泡利矩阵，$J$ 为铁磁耦合强度（$J>0$ 为铁磁，$J<0$ 为反铁磁），$g$ 为横向磁场强度。

通过量子线路演化出状态 $|\psi(t)\rangle$ 后，我们可以在实验中测量以下关键的物理可观测量：

### 1. 平均 Z-磁化强度（Z-Magnetization）
平均自旋磁化代表了铁磁轴方向的经典磁序：
$$\langle M_z \rangle = \frac{1}{N} \sum_{i=1}^{N} \langle \psi(t) | Z_i | \psi(t) \rangle$$
*   **物理意义**：这是体系铁磁相的**序参量（Order Parameter）**。
    *   **强耦合极限（$g \ll J$）**：系统倾向于铁磁自发对称破缺，自旋全部平行（如 $|00...0\rangle$ 或 $|11...1\rangle$），此时 $\langle M_z \rangle \approx 1$。
    *   **强场极限（$g \gg J$）**：横向场占主导，自旋沿 X 方向极化，此时在 Z 方向的磁化强度 $\langle M_z \rangle \approx 0$。
    *   在热力学极限下，系统在 $g/J = 1$ 处发生**量子相变**。

> [!NOTE]
> **在量子计算机中的测量方法**：
> 因为量子计算的计算基底（Computational Basis）就是 Z 的本征态，所以直接对所有比特进行 $Z$-测量即可。
> 假设运行了 $S$ 次实验（shots），第 $i$ 个比特测量得到 $0$ 的次数为 $N_i(0)$，得到 $1$ 的次数为 $N_i(1)$：
> $$\langle Z_i \rangle = \frac{N_i(0) - N_i(1)}{S}$$

---

### 2. 横向 X-极化强度（X-Polarization / Transverse Magnetization）
衡量自旋在横向场方向的响应：
$$\langle M_x \rangle = \frac{1}{N} \sum_{i=1}^{N} \langle \psi(t) | X_i | \psi(t) \rangle$$
*   **物理意义**：代表了由量子横向场引起的量子涨落程度。当 $g \to \infty$ 时，体系处于顺磁相（Paramagnetic Phase），所有自旋被迫指向 X 方向，态为 $|++...+\rangle$，此时 $\langle M_x \rangle \approx 1$。

> [!TIP]
> **在量子计算机中的测量方法**：
> 不能直接测量 $X$，因为硬件只能沿 Z 轴读取。我们必须在测量前对每个比特施加一个 **Hadamard 门（H门）** 进行基底旋转，将 $X$ 轴基底旋转至 $Z$ 轴，然后再进行标准测量。
> 这种测量技术在量子计算中称为**主动基底变换（Active Basis Rotation）**。

---

### 3. 两点自旋关联函数（Two-Point Spin Correlation Function）
两个相距距离为 $r$ 的自旋之间的关联程度：
$$C(r) = \langle Z_i Z_{i+r} \rangle - \langle Z_i \rangle \langle Z_{i+r} \rangle$$
*   **物理意义**：**它在数学上完全等价于量子场论中的格点传播子（Lattice Propagator）！** 
    它代表了量子激发在晶格上的传播和关联。
    *   根据量子统计力学，在临界点 $g/J = 1$ 时，关联长度 $\xi$ 发散，两点关联函数呈幂律衰减 $C(r) \sim r^{-\eta}$（符合共形场论 CFT 的预测）。
    *   在非临界区，关联函数呈指数衰减 $C(r) \sim e^{-r/\xi}$。
    *   在动力学演化中，如果我们创造一个局部激发（比如在一个畴壁域上），$C(r, t)$ 的传播会在时空中形成一个清晰的**“光锥”（Light Cone）**结构，其斜率由 Lieb-Robinson 速度限制。

---

### 4. 纠缠熵（Entanglement Entropy / Von Neumann Entropy）
将体系划分为子系统 $A$（例如左半部分自旋）和子系统 $B$（右半部分），计算子系统 $A$ 的冯·诺伊曼熵：
$$S(A) = -\text{Tr}(\rho_A \ln \rho_A)$$
其中 $\rho_A = \text{Tr}_B(|\psi\rangle\langle\psi|)$ 是子系统 $A$ 的部分密度矩阵。
*   **物理意义**：这是纯粹的量子观测物理量，经典系统绝对不具备。
    *   它定量描述了子系统之间的量子纠缠程度和量子信息的隐形传输。
    *   **纠缠面积律（Area Law）**：对于一维非临界系统，基态的纠缠熵不随系统尺寸增长，而是收敛于常数（即正比于子系统边界大小，1D 的边界就是常数点）。这正是矩阵乘乘积态（MPS）能够高效模拟一维链的物理原因！

---

## 二、 量子计算模拟算法原理

经典计算机无法轻易模拟时间演化算符 $e^{-i H t}$，因为 $H$ 是一个 $2^N \times 2^N$ 的巨大非对角矩阵。量子计算通过以下两个算法实现了高效的动力学演化：

### 1. Trotter-Suzuki 分解（Trotterization）

时间演化算符为 $U(t) = e^{-i H t}$。由于哈密顿量中铁磁相互作用项 $H_Z = -J \sum Z_i Z_{i+1}$ 和横向场项 $H_X = -g \sum X_i$ 互不对易：
$$[H_Z, H_X] \neq 0$$
我们不能直接将指数拆开（即 $e^{-i (H_Z+H_X)t} \neq e^{-i H_Z t} e^{-i H_X t}$）。

**一阶 Trotter 分解公式**将演化时间 $t$ 划分为 $M$ 个极小的时间片 $\Delta t = t/M$：
$$e^{-i H \Delta t} = e^{-i (H_Z + H_X) \Delta t} = e^{-i H_Z \Delta t} e^{-i H_X \Delta t} + \mathcal{O}(\Delta t^2)$$
当 $\Delta t \to 0$（即 $M \to \infty$）时：
$$e^{-i H t} = \lim_{M \to \infty} \left( e^{-i H_Z \Delta t} e^{-i H_X \Delta t} \right)^M$$

#### 在量子线路中的逻辑门映射：
我们将演化算符拆解为可以在量子芯片上直接执行的局部量子门：
1.  **单体项演化（横向场）**：
    $$e^{-i g \Delta t X_i} = R_x(2g\Delta t)_i$$
    这直接对应于每个比特上的旋转门 **$R_x$**。
2.  **双体项演化（自旋耦合）**：
    $$e^{i J \Delta t Z_i Z_{i+1}} = R_{zz}(-2J\Delta t)_{i, i+1}$$
    这个两比特旋转门通过标准量子门分解为：一个 $CNOT$ 门，在目标位施加 $R_z(-2J\Delta t)$，然后再接一个 $CNOT$ 门。
    
    ```text
    q_i:   ───■────────────────■───
            │  ┌────────────┐  │   
    q_i+1: ───X──┤ Rz(θ_Ising) ├─X───
    ```
    *(其中 $\theta_{Ising} = -2J\Delta t$)*

---

### 2. MPS（矩阵乘积态）张量网络模拟原理

对于大体系模拟（如您运行的 $N=60$ spins），Qiskit Aer 并没有维护 $2^{60}$ 维度的态矢量，而是使用 **MPS（Matrix Product State）** 进行高效压缩。

#### ① 数学形式
MPS 将状态矢量的每一个振幅系数 $c_{s_1 s_2 ... s_N}$ 分解为一连串局部小矩阵的乘积：
$$| \psi \rangle = \sum_{s_1, ..., s_N} \text{Tr}(A_1^{s_1} A_2^{s_2} \cdots A_N^{s_N}) | s_1 s_2 \cdots s_N \rangle$$
*   其中 $s_i \in \{0, 1\}$ 是第 $i$ 个物理量子比特的基态。
*   $A_i^{s_i}$ 是一个 $\chi \times \chi$ 的复数矩阵。$\chi$ 被称为**键维（Bond Dimension）**或 Schmidt 秩。

#### ② 局域幺正更新（Local Gate Application）
当我们在第 $i$ 和 $i+1$ 个量子比特上施加一个 CNOT 门或 Rzz 门时，MPS 不需要全局更新：
1.  它仅仅将局域的两个张量 $A_i$ 和 $A_{i+1}$ 进行收缩合并，形成一个较大的临时张量。
2.  在这个合并张量上进行 **奇异值分解（SVD）**。
3.  通过保留前 $\chi$ 个最大的奇异值，丢弃极小的值（截断），重新将其分裂回两个新的张量 $A'_i$ 和 $A'_{i+1}$。

#### ③ 为什么在大系统下 MPS 能够保证高精度？
*   一维量子链在局域 Trotter 步演化下，**纠缠只能以有限的 Lieb-Robinson 速度从接触点向外传播**。
*   因此，Schmidt 谱（奇异值 $\lambda_i$）随着索引呈指数衰减。
*   只要我们截断时保留的键维 $\chi$ 足够大（例如 $\chi=128$ 或 $256$），丢弃的奇异值之和（即**截断误差**）就会极其微小（$< 10^{-8}$）。
*   这使得 MPS 能够在保证**近乎物理精确**的同时，将计算复杂度从状态矢量的 $\mathcal{O}(2^N)$ 彻底降低到张量网络的 $\mathcal{O}(N \cdot \chi^3)$！

---

## 总结：QFT 与 TFIM 之间的深层对偶性

既然您正在探究量子场论，TFIM 实际上可以在一维空间中通过 **Jordan-Wigner 变换** 精确映射为一维的**无质量/有质量狄拉克费米子场**！
*   铁磁相和顺磁相之间的量子临界点，在场论中正好对应着一个**无质量的自旋-1/2  Majorana 费米子共形场论（CFT）**。
*   您手写的谐振子对角化（`05_anharmonic_osc.py`）以及格点传播子计算（`07_phi4_lattice_propagator.py`），与我们用量子线路模拟 TFIM 并在大系统上用张量网络收缩获取自旋关联函数 $C(r, t)$，在物理内核上是**完全相通的**。量子计算提供了一种全新的、在指数级希尔伯特空间中操纵这些场算符和波函数的终极手段。
