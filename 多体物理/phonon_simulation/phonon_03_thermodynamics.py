"""
Phonon Theory Script 03: Debye and Einstein Models for Solid Specific Heat & Thermodynamics
========================================================================================
Formula to Code Mapping:
- Eq (3.1) Bose-Einstein Phonon Occupation Number:
  <n(\omega)> = \frac{1}{e^{\hbar\omega / k_B T} - 1}
  -> Function `bose_einstein_occupation(omega, T)`
- Eq (3.2) Einstein Specific Heat:
  C_V^E(T) = 3 N k_B (\Theta_E / T)^2 * \frac{e^{\Theta_E/T}}{(e^{\Theta_E/T} - 1)^2}
  -> Function `c_v_einstein(T, Theta_E)`
- Eq (3.3) Debye Specific Heat:
  C_V^D(T) = 9 N k_B (T / \Theta_D)^3 \int_0^{\Theta_D/T} \frac{x^4 e^x}{(e^x - 1)^2} dx
  -> Function `c_v_debye(T, Theta_D)`
- Eq (3.4) Low-T Debye Limit:
  C_V^D(T) \approx \frac{12 \pi^4}{5} N k_B (T / \Theta_D)^3
  -> Function `c_v_debye_low_T(T, Theta_D)`
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150

kB = 1.0  # Normalized units (kB = 1)
N = 1.0   # Normalized 1 mol / 1 atom unit

def bose_einstein_occupation(omega, T):
    """ Formula (3.1): Bose-Einstein Distribution Function """
    if T <= 0:
        return 0.0
    x = omega / T
    # Avoid overflow in exp(x)
    x = np.clip(x, 1e-10, 500)
    return 1.0 / (np.exp(x) - 1.0)

def c_v_einstein(T, Theta_E=1.0):
    """ Formula (3.2): Einstein Specific Heat """
    T = np.atleast_1d(T)
    Cv = np.zeros_like(T, dtype=float)
    for i, t in enumerate(T):
        if t <= 0:
            Cv[i] = 0.0
        else:
            x = Theta_E / t
            if x > 200:
                Cv[i] = 0.0
            else:
                exp_x = np.exp(x)
                Cv[i] = 3.0 * N * kB * (x**2) * exp_x / ((exp_x - 1.0)**2)
    return Cv

def debye_integrand(x):
    """ Integrand for Debye model: x^4 * e^x / (e^x - 1)^2 """
    if x < 1e-6:
        return x**2  # Series expansion for small x to prevent 0/0
    if x > 200:
        return 0.0
    exp_x = np.exp(x)
    return (x**4 * exp_x) / ((exp_x - 1.0)**2)

def c_v_debye(T, Theta_D=1.0):
    """ Formula (3.3): Debye Specific Heat Integrator """
    T = np.atleast_1d(T)
    Cv = np.zeros_like(T, dtype=float)
    for i, t in enumerate(T):
        if t <= 0:
            Cv[i] = 0.0
        else:
            x_max = Theta_D / t
            integral, _ = quad(debye_integrand, 0, x_max)
            Cv[i] = 9.0 * N * kB * ((t / Theta_D)**3) * integral
    return Cv

def c_v_debye_low_T(T, Theta_D=1.0):
    """ Formula (3.4): T^3 Low-Temperature Asymptote """
    return (12.0 * np.pi**4 / 5.0) * N * kB * ((T / Theta_D)**3)

def plot_phonon_thermodynamics():
    Theta = 1.0 # Temperature normalized to characteristic scale
    T_range = np.linspace(0.001, 2.5, 300)

    Cv_E = c_v_einstein(T_range, Theta_E=Theta)
    Cv_D = c_v_debye(T_range, Theta_D=Theta)
    Cv_D_low = c_v_debye_low_T(T_range, Theta_D=Theta)

    fig = plt.figure(figsize=(14, 6))

    # --- Panel 1: Specific Heat Comparison C_V(T) ---
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(T_range / Theta, Cv_D, 'b-', lw=2.5, label='德拜模型 $C_V^{Debye}(T)$')
    ax1.plot(T_range / Theta, Cv_E, 'r--', lw=2.5, label='爱因斯坦模型 $C_V^{Einstein}(T)$')
    ax1.axhline(3.0 * N * kB, color='black', linestyle=':', label='杜隆-珀蒂经典极限 ($3Nk_B$)')
    
    # Low-T zoom inset or curve
    ax1.plot(T_range / Theta, Cv_D_low, 'g-.', lw=1.8, label=r'低温 $T^3$ 渐近律 $\frac{12\pi^4}{5}Nk_B(T/\Theta_D)^3$')

    ax1.set_xlabel('温度 $T / \Theta$', fontsize=12)
    ax1.set_ylabel('等容比热 $C_V / (N k_B)$', fontsize=12)
    ax1.set_title('图 3a: 固体声子比热容模型对比 (Debye vs Einstein)', fontsize=14, pad=12)
    ax1.set_xlim([0, 2.5])
    ax1.set_ylim([0, 3.5])
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend(fontsize=10, loc='lower right')

    # --- Panel 2: Low Temperature T^3 Behavior Log-Log Plot ---
    ax2 = fig.add_subplot(1, 2, 2)
    T_low = np.logspace(-2, -0.7, 100) # T/Theta from 0.01 to 0.2
    Cv_D_log = c_v_debye(T_low, Theta_D=Theta)
    Cv_E_log = c_v_einstein(T_low, Theta_E=Theta)
    Cv_low_log = c_v_debye_low_T(T_low, Theta_D=Theta)

    ax2.loglog(T_low, Cv_D_log, 'b-', lw=2.5, label='德拜模型 $C_V^{Debye}$')
    ax2.loglog(T_low, Cv_E_log, 'r--', lw=2.5, label='爱因斯坦模型 $C_V^{Einstein}$')
    ax2.loglog(T_low, Cv_low_log, 'g-.', lw=2, label=r'严格 $T^3$ 幂律斜率')

    ax2.set_xlabel('温度 $\ln(T / \Theta_D)$', fontsize=12)
    ax2.set_ylabel('等容比热 $\ln(C_V)$', fontsize=12)
    ax2.set_title('图 3b: 低温下德拜 $T^3$ 律与爱因斯坦指数冻结的双对数图', fontsize=14, pad=12)
    ax2.grid(True, which="both", linestyle=':', alpha=0.6)
    ax2.legend(fontsize=10, loc='upper left')

    # Annotation for freezing
    ax2.annotate('爱因斯坦模型:\n单频光子模低温呈指数冻结 $e^{-\Theta_E/T}$',
                 xy=(0.04, c_v_einstein(0.04, 1.0)[0]),
                 xytext=(0.02, 1e-4),
                 arrowprops=dict(facecolor='red', shrink=0.08, width=1, headwidth=5),
                 fontsize=9, bbox=dict(boxstyle="round,pad=0.3", fc="MistyRose", ec="r", lw=1))

    ax2.annotate('德拜模型:\n低频声学支连续激发 $\propto T^3$',
                 xy=(0.05, c_v_debye(0.05, 1.0)[0]),
                 xytext=(0.08, 1e-2),
                 arrowprops=dict(facecolor='blue', shrink=0.08, width=1, headwidth=5),
                 fontsize=9, bbox=dict(boxstyle="round,pad=0.3", fc="LightCyan", ec="b", lw=1))

    plt.tight_layout()
    plt.savefig('phonon_03_thermodynamics.png', dpi=300)
    print("Saved figure: phonon_03_thermodynamics.png")

if __name__ == '__main__':
    plot_phonon_thermodynamics()
