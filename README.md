# qft_coding — Learning QFT Through Computation

从谐振子到 Feynman 图，每个概念都有对应的数值验证。

**路线：** QHO → 耦合谐振子链（自由场）→ 非谐振子（微扰论）→ 格点 φ⁴（tadpole）→ ...

**文档：** `qft_basics.tex`（QFT 基础，简体中文，含完整推导和数值对照）

详细路线图见 [ROADMAP.md](ROADMAP.md)。

---

## 练习列表

### 基础 QFT（从 QM 出发）

| File | Topic | Status |
|------|-------|--------|
| `04_coupled_qho.py` | 耦合谐振子链 → 色散关系 ω(k) = √(m² + 4g sin²(k/2)) | ✅ |
| `05_anharmonic_osc.py` | 单非谐振子：精确对角化 vs Rayleigh-Schrödinger 微扰 | ✅ |
| `06_phi4_lattice_tadpole.py` | 格点 φ⁴：tadpole 模式求和 → Feynman 图 | ✅ |

### 有限温专题

| File | Topic | Status |
|------|-------|--------|
| `01_tadpole_thermal_mass.py` | Tadpole & thermal mass m²_th = λT²/24 | ✅ |
| `02_spectral_function.py` | Spectral function & thermal broadening | ✅ |
| `03_HTL_plasmon.py` | HTL self-energy & plasmon dispersion | ✅ |

---

## 使用

```bash
python 06_phi4_lattice_tadpole.py    # 运行格点 φ⁴ 数值验证
```

依赖：numpy, scipy, matplotlib

## 参考

- Peskin & Schroeder, *An Introduction to Quantum Field Theory*
- Kapusta & Gale, *Finite-Temperature Field Theory*
