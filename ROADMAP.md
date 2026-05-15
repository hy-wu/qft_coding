# qft_coding — Learning QFT Through Computation

从量子力学到场论，从谐振子到 Feynman 图。
每个概念都有对应的 Python 数值验证和 TeX 文档。

---

## Roadmap

### Phase 0：从谐振子到算符场 (QHO → Field Operator)

**目标：** 搞清楚算符场 $\hat\phi(x)$ 到底是什么。

| Step | Topic | Code | Status |
|------|-------|------|--------|
| 0.0 | 单量子谐振子：产生湮灭算符、Fock 空间、零点能 | — (conceptual) | ✅ |
| 0.1 | 耦合谐振子链 → 简正模式 → 色散关系 $\omega_k^2 = m^2 + 4g\sin^2(k/2)$ | `04_coupled_qho.py` | ✅ |
| 0.2 | 连续极限 → 场算符 $\hat\phi(x)$，传播子的起源 | `qft_basics.tex` | ✅ |

**核心结论：** 每个 $k$ 模式是独立谐振子。粒子 = 场的量子化激发。
真空 = 所有模式处于基态。真空涨落 → 传播子。

---

### Phase 1：自由标量场 (Free Scalar Field)

| Step | Topic | Code | Status |
|------|-------|------|--------|
| 1.1 | 耦合 QHO 链对角化，数值色散关系 | `04_coupled_qho.py` | ✅ |

---

### Phase 2：相互作用与微扰论 (Interactions & PT)

**目标：** 从 $\lambda x^4$ 到 $\phi^4$ tadpole 图。

| Step | Topic | Code | Status |
|------|-------|------|--------|
| 2.0 | 单非谐振子：精确对角化 vs 1阶/2阶 Rayleigh-Schrödinger 微扰 | `05_anharmonic_osc.py` | ✅ |
| 2.1 | 格点 $\phi^4$：tadpole 的模式求和 → Feynman 图 | `06_phi4_lattice_tadpole.py` | ✅ |
| 2.2 | 格点传播子：等时两点关联函数 $C_{ij} = \langle 0\vert\phi_i\phi_j\vert 0\rangle$ 与质量重整化 | `07_phi4_lattice_propagator.py` | ✅ |
| 2.3 | 动量空间传播子与自能 $\Sigma(k)$: Dyson 方程, tadpole 自能 | `08_phi4_self_energy.py` | ✅ |

**Phase 2 核心结果：**
- $\langle 0|x^4|0\rangle = 3/(4\omega^2)$ 是单个谐振子的 tadpole
- 格点 $\phi^4$ 中 $\langle 0|\phi_i^4|0\rangle = 3(\frac1N\sum_k\frac1{2\omega_k})^2$
- 解析公式与数值精确对角化一致（误差 $2.6\times10^{-18}$）

---

### Phase 3：散射与截面 (Scattering & Cross Sections)

| Step | Topic | Code | Status |
|------|-------|------|--------|
| 3.0 | LSZ 约化公式与 $S$ 矩阵 | — | 📋 planned |
| 3.1 | $\phi^4$ 的 $2\to2$ 散射：树图截面 | — | 📋 planned |
| 3.2 | 单圈修正：$s$, $t$, $u$ 通道 | — | 📋 planned |

---

### Finite-T 专题 (早期练习)

| # | Topic | Code | Status |
|---|-------|------|--------|
| 01 | Tadpole & thermal mass: $m_{\rm th}^2 = \lambda T^2/24$ | `01_tadpole_thermal_mass.py` | ✅ |
| 02 | Spectral function: Breit-Wigner, thermal width $\Gamma \sim \lambda^2 T$ | `02_spectral_function.py` | ✅ |
| 03 | HTL resummation: plasmon dispersion, Debye mass | `03_HTL_plasmon.py` | ✅ |

---

## 文档

- **`qft_basics.tex`** — QFT 基础：从谐振子到场算符、微扰论、格点 $\phi^4$（14 页，简体中文）
- **`paper.tex`** — 有限温场论：tadpole、谱函数、HTL（简体中文）

---

## 依赖

- numpy, scipy, matplotlib
- xelatex（编译 TeX → PDF）

---

## 参考

- Peskin & Schroeder, *An Introduction to Quantum Field Theory*
- Kapusta & Gale, *Finite-Temperature Field Theory*
- Laine & Vuorinen, *Basics of Thermal Field Theory*
