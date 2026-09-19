# 微扰展开的展开参数与有效性 (What Is Expanded, and Why It Works)

> **主要目标**：回答两个问题 —— Dyson 级数展开的到底是**耦合**还是 **$\hbar$**？
> 以及**为什么这个展开有效**？
>
> 简短答案：展开的是**耦合**；$\hbar$ 数的是**圈数**；在**固定外腿数**下两者只差一个平移。
> 有效性来自"零阶可解 + 渐近级数 + 耦合小"，而**不是**来自任何收敛性 ——
> 收敛性根本不存在。
>
> 配套阅读：[`wick_expansion.md`](wick_expansion.md)、[`feynman_diagrams.md`](feynman_diagrams.md)、
> [`qft_before_computers.md`](qft_before_computers.md)。

---

## 一、 展开参数是耦合，不是 $\hbar$

Dyson 级数第 $n$ 项的系数是

$$
\frac{(-i)^n}{n!}\Big(\frac{\lambda}{4!}\Big)^n
$$

幂次是 $\lambda^n$ —— 展开参数是**耦合常数**，与 $\hbar$ 无关（全文取 $\hbar=1$）。

QM 里同理：$U_I(t,t_0)=T\exp\big[-\frac{i}{\hbar}\int H_I\,dt\big]$，展开的幂次是 $H_I$ 的强度，
即"相互作用能 / 能级间距" $\ll 1$。**也不是 $\hbar$ 展开。**

---

## 二、 把 $\hbar$ 恢复：它数的是圈数

在路径积分里恢复 $\hbar$：

$$
Z=\int\mathcal D\phi\;\exp\Big[\frac{i}{\hbar}S[\phi]\Big],
\qquad
S=\int d^4x\Big[\tfrac12(\partial\phi)^2-\tfrac12m^2\phi^2-\frac{\lambda}{4!}\phi^4\Big]
$$

逐项数 $\hbar$ 的幂：

| 元素 | 来源 | $\hbar$ 幂次 |
|---|---|---|
| 传播子 | 二次项的高斯积分 $\to K^{-1}$ | $+1$ |
| 顶点 | 展开 $e^{\frac{i}{\hbar}S_{\rm int}}$ | $-1$ |
| 外线 | 渐近态归一化 | $0$ |

一张有 $V$ 个顶点、$I$ 条内线的图：

$$
\hbar^{\,I}\cdot\hbar^{-V}=\hbar^{\,I-V}\;\overset{L=I-V+1}{=}\;\hbar^{\,L-1}
$$

> **$\hbar$ 数的是圈数，不是顶点数。每多一个圈，多一个 $\hbar$。**

---

## 三、 固定外腿时，两者只差一个平移

对 $\phi^4$，数"线头"：

$$
4V=\underbrace{E}_{\text{外腿各 1 个}}+\underbrace{2I}_{\text{内线两端}},
\qquad L=I-V+1
$$

联立解出

$$
\boxed{\;V=\frac{E}{2}+L-1\;}\qquad\Longleftrightarrow\qquad L=V-\frac{E}{2}+1
$$

固定 $E$ 时 $V$ 与 $L$ **一一对应**：

| 过程 | $E$ | 树图 $L{=}0$ | 1 圈 $L{=}1$ | 2 圈 $L{=}2$ |
|---|---|---|---|---|
| 自能 | 2 | $V=0$（自由传播子） | $V=1$（tadpole） | $V=2$ |
| $2\to2$ | 4 | $V=1$（接触项） | $V=2$ | $V=3$ |
| 真空 | 0 | — | — | $V=1$（八字形） |

所以对 $2\to2$ 散射：

$$
\underbrace{\lambda^1}_{\text{Dyson 第 1 阶}}\longleftrightarrow\underbrace{\hbar^{0}}_{\text{树图}},\qquad
\underbrace{\lambda^2}_{\text{第 2 阶}}\longleftrightarrow\underbrace{\hbar^{1}}_{\text{1 圈}},\qquad
\underbrace{\lambda^3}_{\text{第 3 阶}}\longleftrightarrow\underbrace{\hbar^{2}}_{\text{2 圈}}
$$

**同一件事，两种数法。** 但注意这个对应**只在固定外腿数时成立**。若把不同 $E$ 的过程混在一起
比较"阶数"（例如拿 $2\to2$ 树图与自能 tadpole 比），对应就断了。

---

## 四、 真正"按 $\hbar$ 展开"的是圈展开 / 半经典展开

它来自路径积分的**驻相法**（stationary phase）：

$$
Z=\int\mathcal D\phi\;e^{\frac{i}{\hbar}S[\phi]}
\;\xrightarrow{\ \hbar\to0\ }\;
\underbrace{e^{\frac{i}{\hbar}S_{\rm cl}}}_{\text{经典解}}\times\underbrace{\Big(\text{Gaussian 涨落}\Big)}_{\text{一圈}}\times\cdots
$$

有效作用量的展开正是这个：

$$
\Gamma[\phi]=\underbrace{S_{\rm cl}[\phi]}_{\hbar^{-1}}
+\underbrace{\hbar\,\Gamma_{1\rm loop}[\phi]}_{L=1}
+\underbrace{\hbar^2\,\Gamma_{2\rm loop}[\phi]}_{L=2}+\cdots
$$

**QM 与 QFT 的平行结构**：

| | 耦合展开 | $\hbar$ 展开 |
|---|---|---|
| QM | Dyson / Rayleigh–Schrödinger 微扰（$\lambda$） | **WKB / 半经典** |
| QFT | Dyson 级数（顶点数 $V$） | **圈展开**（圈数 $L$） |

两个是**不同的**展开；只是在 QFT 里、固定外腿时被第三节的关系绑在了一起。

---

## 五、 为什么有效

### 5.1 为什么"能用"

因为**零阶可解**：自由场是精确可解的高斯理论，其 Hilbert 空间（Fock 空间）就是我们熟悉的
"粒子"图像。相互作用只是在这个可解图像上加一点扰动，所以按扰动强度展开，前几项就能
抓住主要物理。

### 5.2 但它**不收敛** —— 是渐近级数

**Dyson 反证（1952）**：假设 $F(\lambda)=\sum_n a_n\lambda^n$ 在某个 $\lambda_0>0$ 收敛。
则它在 $|\lambda|<\lambda_0$ 内解析，因此对**小负耦合**也解析。但 $\lambda<0$ 时势能

$$
V(x)=\tfrac12m^2x^2+\frac{\lambda}{4!}x^4\qquad(\lambda<0)
$$

**无下界** —— 没有稳定基态，$E_0(\lambda)$ 根本不是解析函数。矛盾。**故收敛半径为零。**

定量地：第 $n$ 阶图的个数随 $n$ 按 $n!$ 增长，于是

$$
a_n\;\sim\;C^{\,n}\,n!\qquad\Longrightarrow\qquad R_{\rm conv}=0
$$

### 5.3 那它靠什么"有效"：最优截断

虽然不收敛，但它是**渐近级数**：截断到最优阶 $n^*$ 时误差最小，且

$$
n^*\sim\frac{1}{\lambda},
\qquad
\text{最小误差}\;\sim\;e^{-c/\lambda}
$$

$e^{-c/\lambda}$ 是**非微扰**效应（瞬子 / 隧穿）。所以微扰论的精度上限由非微扰效应决定 ——
有些物理量（QCD 的禁闭、拓扑荷）微扰论永远给不出来。

### 5.4 有效性的实际判据：无量纲耦合要小

| 理论 | 耦合 | 微扰论 |
|---|---|---|
| QED（低能） | $\alpha=1/137$ | 极准（$a_e$ 十位有效数字） |
| QCD（高能） | $\alpha_s\ll1$（渐近自由） | 准 |
| QCD（$\lesssim1$ GeV） | $\alpha_s\sim1$ | **失效** → 格点 |
| 4 维 $\phi^4$ | $\lambda$ 跑动，UV 有 Landau 极点 | IR 弱耦合可用；UV 不干净 |

---

## 六、 两个更深的坑（不影响使用，但值得知道）

1. **Renormalon**：即使把 $n!$ 发散按 Borel 求和，4 维 $\phi^4$ 的 IR renormalon 也会带来
   **同号**阶乘增长，使 Borel 和本身不唯一。严格的理论定义只能靠**非微扰**方法。
2. **Haag 定理**：相互作用绘景在 QFT 里**严格不存在** —— 自由理论与相互作用理论不幺正等价。
   所以 Dyson 级数从一开始就是个**形式展开**。现代观点：先用格点 / 构造性方法**定义**理论，
   微扰论只是它的渐近展开。

---

## 七、 与本仓库的对应

这不是抽象讨论 —— Phase 0–2 就是这条逻辑的实证：

- **`05_anharmonic_osc.py`**：ROADMAP 里写着"小 $\lambda$（$<0.1$）时 1 阶和 2 阶微扰都很好；
  **大 $\lambda$ 时微扰论发散**，精确对角化仍然是可靠的"。这正是渐近级数的行为 ——
  小耦合好、大耦合坏，而且坏得无法靠加阶数挽救。
- **`06_phi4_lattice_tadpole.py` / `07_phi4_lattice_propagator.py`**：格点精确对角化就是那个
  **非微扰定义**。微扰论是它的渐近展开，而不是反过来。
- 所以 Phase 2 用格点去"验证"微扰论的方向是对的：**微扰论本身没有独立意义**，
  它的合法性来自与格点（或实验）的一致。

---

## 一句话总结

Dyson 级数展开的是**耦合**；$\hbar$ 数的是**圈**；固定外腿时两者相差一个平移，
所以"耦合展开 $=$ 圈展开 $=$ $\hbar$ 展开"在那个意义上成立。
它的有效性来自**零阶可解 + 渐近级数 + 耦合小**，而非任何收敛性。
