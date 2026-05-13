"""
有限溫度 φ⁴ 理論 — 單圈 tadpole 與熱質量
=========================================
場論計算 → 數值實現 → 物理結果

流程:
  1. 寫出松原 propagator: Δ(iωₙ, k) = 1/(ωₙ² + k² + m₀²)
  2. 單圈 tadpole 自能: Π = λ/2 T Σₙ ∫ d³k/(2π)³ Δ(iωₙ, k)
  3. 解析做松原和 → 數值做動量積分
  4. 畫出 m_th²(T) 並和高溫極限 λT²/24 比較

背景:
  QGP 中的準粒子在高溫下獲得熱質量，這是德拜屏蔽的場論表現。
  對於 φ⁴ 理論，高溫極限下 m_th² = λT²/24。
  QCD 中膠子的熱質量 m_D² = (g²T²/3)(N_c + N_f/2)，是等離子體頻率的來源。
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

# ================================================================
# 物理參數 (自然單位 ℏ = c = k_B = 1)
# ================================================================
LAMBDA = 1.0     # φ⁴ 耦合常數
M0 = 0.1         # 零溫質量 (紅外調節)

# ================================================================
# 第一步: 解析松原和 → 數值積分
#
# 松原頻率: ωₙ = 2πnT (玻色子)
# Tadpole 自能:
#
#   Π = λ/2 * T Σ_{n=-∞}^{∞} ∫ d³k/(2π)³  1/(ωₙ² + Eₖ²)
#     其中 Eₖ = √(k² + m₀²)
#
# 對 n 求和是標準的 (可在 Kapusta & Gale, Laine & Vuorinen 中找到):
#
#   T Σₙ 1/(ωₙ² + E²) = (1/2E)[1 + 2n_B(E)]
#                        ─┬─────  ──┬───
#                      T=0 真空     熱部分
#
# 零溫部分是發散的，但在計算熱質量時，我們減去 T=0 部分
# (重整化後只保留溫度貢獻):
#
#   Π_th(T) = λ/2 ∫ d³k/(2π)³  n_B(Eₖ)/Eₖ
#
# 球坐標: d³k = 4π k² dk
#
#   Π_th(T) = λ/(4π²) ∫₀^∞ dk  k² n_B(Eₖ)/Eₖ
# ================================================================

def bose_einstein(E, T):
    """玻色-愛因斯坦分佈 n_B(E) = 1/(e^{E/T} - 1)"""
    if T <= 0:
        return 0.0
    # 防止 exp 溢出
    if E / T > 100:
        return 0.0
    return 1.0 / (np.exp(E / T) - 1.0)


def tadpole_integrand(k, T):
    """被積函數: 4π k² n_B(Eₖ)/Eₖ  (球坐標因子已包含)"""
    Ek = np.sqrt(k**2 + M0**2)
    return 4.0 * np.pi * k**2 * bose_einstein(Ek, T) / Ek


def thermal_self_energy(T):
    """
    數值計算 tadpole 對自能的溫度貢獻。

    積分截斷: 動量 ~ 幾個 T 的熱區貢獻最大，
              取 cutoff = max(20T, 10) 足夠。
    """
    if T <= 1e-10:
        return 0.0
    cutoff = max(20.0 * T, 10.0)
    result, err = integrate.quad(tadpole_integrand, 0, cutoff,
                                 args=(T,), limit=200, epsabs=1e-8, epsrel=1e-6)
    return LAMBDA / 2.0 * result / (2 * np.pi)**3


def thermal_mass_sq_HT(T):
    """高溫極限: m_th² = λT²/24"""
    return LAMBDA * T**2 / 24.0


# ================================================================
# 第二步: 掃描溫度並計算
# ================================================================
# T 從 0.1 m₀ 到 ~30 m₀，涵蓋低溫到高溫區
temperatures = np.logspace(np.log10(0.3), np.log10(30.0), 16)
m2_numerical = []
m2_HT_limit = []

print(f"{'T/m0':>8}  {'m_th^2(num)':>14}  {'m_th^2(HTL)':>14}  {'ratio':>8}")
print("-" * 50)

for T in temperatures:
    m2_num = thermal_self_energy(T)
    m2_ht = thermal_mass_sq_HT(T)
    m2_numerical.append(m2_num)
    m2_HT_limit.append(m2_ht)
    ratio = m2_num / m2_ht if m2_ht > 0 else 0
    print(f"{T/M0:>8.2f}  {m2_num:>14.6f}  {m2_ht:>14.6f}  {ratio:>8.4f}")

# ================================================================
# 第三步: 畫圖
# ================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# 左圖: 熱質量 vs 溫度 (對數坐標)
ax1.plot(temperatures / M0, m2_numerical, 'o-', color='#1f77b4',
         label='數值積分', linewidth=2, markersize=6)
ax1.plot(temperatures / M0, m2_HT_limit, '--', color='#d62728',
         label='高溫極限 λT²/24', linewidth=2)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel(r'$T / m_0$', fontsize=14)
ax1.set_ylabel(r'$m_{\rm th}^2(T)$', fontsize=14)
ax1.set_title(r'$\phi^4$ Thermal Mass vs Temperature', fontsize=13)
ax1.legend(fontsize=12)
ax1.grid(True, alpha=0.3)

# 右圖: 與高溫極限的比值
ratio_arr = np.array(m2_numerical) / np.array(m2_HT_limit)
ax2.semilogx(temperatures / M0, ratio_arr, 's-', color='#2ca02c',
             linewidth=2, markersize=6)
ax2.axhline(y=1.0, color='gray', linestyle=':', alpha=0.7)
ax2.set_xlabel(r'$T / m_0$', fontsize=14)
ax2.set_ylabel(r'$m_{\rm th}^2 / (\lambda T^2/24)$', fontsize=14)
ax2.set_title(r'Deviation from High-T Limit', fontsize=13)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('01_thermal_mass.png', dpi=150)
print(f"\nSaved: 01_thermal_mass.png")

# ================================================================
# 第四步: 驗證 — 高溫下 m_th²/T² 趨於常數
# ================================================================
print("\n--- High-T limit verification ---")
for T in temperatures[-4:]:
    print(f"T = {T:.2f}  m_th^2/T^2 = {thermal_self_energy(T) / T**2:.6f}  "
          f"(theory: lambda/24 = {LAMBDA/24:.6f})")
