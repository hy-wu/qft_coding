r"""
Phonon Theory Script 04: Quantum Field Fluctuation and Second Quantization Observables
=====================================================================================
Formula to Code Mapping:
- Eq (4.1) Second Quantized Field Operator:
  \hat{u}_n = \sum_k \sqrt{\frac{\hbar}{2 N M \omega_k}} ( a_k e^{ikna} + a_k^\dagger e^{-ikna} )
- Eq (4.2) Mean Square Displacement Fluctuation:
  < \hat{u}_n^2 > = \frac{\hbar}{2 N M} \sum_{k} \frac{1}{\omega_k} \coth\left( \frac{\hbar \omega_k}{2 k_B T} \right)
  -> Function `mean_square_displacement(T, N, M, C, a)`
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150

def compute_quantum_fluctuations(N=200, T_list=np.linspace(0, 2.0, 100), hbar=1.0, kB=1.0, M=1.0, C=1.0, a=1.0):
    """
    Formula (4.2): Evaluates quantum thermal mean-square displacement <u_n^2>
    as a function of temperature T.
    """
    n_k = np.arange(1, N // 2)
    k_vals = 2.0 * np.pi * n_k / (N * a)
    omega_k = 2.0 * np.sqrt(C / M) * np.sin(k_vals * a / 2.0)

    u2_mean = []
    u2_zero_point = 0.5 * hbar / (N * M) * np.sum(1.0 / omega_k)

    for T in T_list:
        if T == 0:
            coth_factor = np.ones_like(omega_k)
        else:
            x = hbar * omega_k / (2.0 * kB * T)
            x_clipped = np.clip(x, 1e-6, 100)
            coth_factor = np.cosh(x_clipped) / np.sinh(x_clipped)
            
        u2 = (hbar / (N * M)) * np.sum((1.0 / omega_k) * coth_factor)
        u2_mean.append(u2)

    return np.array(T_list), np.array(u2_mean), u2_zero_point, k_vals, omega_k

def plot_quantum_fluctuations():
    N_atoms = 300
    T_range = np.linspace(0, 2.0, 100)
    T_eval, u2_eval, u2_zp, k_vals, omega_k = compute_quantum_fluctuations(N=N_atoms, T_list=T_range)

    fig = plt.figure(figsize=(14, 6))

    # --- Panel 1: Mean Square Quantum & Thermal Fluctuation vs Temperature ---
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(T_eval, u2_eval, 'b-', lw=2.5, label=r'总位移均方涨落 $\langle \hat{u}_n^2 \rangle(T)$')
    ax1.axhline(u2_zp, color='crimson', linestyle='--', lw=2, label=r'基态零点能量子涨落 $\langle \hat{u}_n^2 \rangle_{T=0}$')

    ax1.annotate(f'量子零点涨落 (Zero-Point Fluctuation)\n' + r'$\langle \hat{u}_n^2 \rangle_{T=0} = $' + f'{u2_zp:.4f}',
                 xy=(0, u2_zp), xytext=(0.25, u2_zp * 1.5),
                 arrowprops=dict(facecolor='red', shrink=0.08, width=1.5, headwidth=6),
                 fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="MistyRose", ec="r", lw=1))

    ax1.plot(T_eval[T_eval > 0.5], T_eval[T_eval > 0.5] * 0.45 + u2_zp, 'k:', label=r'高温经典能量均分线性渐近线')

    ax1.set_xlabel(r'温度 $T / (\hbar \omega_{max} / k_B)$', fontsize=12)
    ax1.set_ylabel(r'原子位移均方涨落 $\langle \hat{u}_n^2 \rangle$', fontsize=12)
    ax1.set_title(r'图 4a: 二次量子化场算符的量子零点涨落与热涨落', fontsize=14, pad=12)
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend(fontsize=10, loc='upper left')

    # --- Panel 2: Simulated Real-space Quantum Displacement Waveform Snapshot ---
    ax2 = fig.add_subplot(1, 2, 2)
    np.random.seed(42)
    n_sites = np.arange(100)

    phase_zp = np.random.uniform(0, 2*np.pi, len(k_vals))
    amp_zp = np.sqrt(1.0 / omega_k)
    u_snapshot_zp = np.zeros(len(n_sites))
    
    phase_th = np.random.uniform(0, 2*np.pi, len(k_vals))
    coth_th = np.cosh(omega_k / 2.0) / np.sinh(omega_k / 2.0)
    amp_th = np.sqrt((1.0 / omega_k) * coth_th)
    u_snapshot_th = np.zeros(len(n_sites))

    for ik, k in enumerate(k_vals):
        u_snapshot_zp += amp_zp[ik] * np.cos(k * n_sites + phase_zp[ik])
        u_snapshot_th += amp_th[ik] * np.cos(k * n_sites + phase_th[ik])

    u_snapshot_zp /= np.std(u_snapshot_zp)
    u_snapshot_th /= np.std(u_snapshot_th)

    ax2.plot(n_sites, u_snapshot_zp, 'r.-', alpha=0.8, label=r'基态量子零点波动快照 ($T=0$)')
    ax2.plot(n_sites, u_snapshot_th + 3.5, 'b.-', alpha=0.8, label=r'热激发态波动快照 ($T > 0$)')

    ax2.set_xlabel(r'原子编号 $n$', fontsize=12)
    ax2.set_ylabel(r'标准化位移偏离 $u_n$', fontsize=12)
    ax2.set_title(r'图 4b: 实空间原子位移算符 $\hat{u}_n(t)$ 的量子场快照', fontsize=14, pad=12)
    ax2.set_yticks([0, 3.5])
    ax2.set_yticklabels(['基态 (T=0)', '热激发态'])
    ax2.grid(True, linestyle=':', alpha=0.5)
    ax2.legend(fontsize=10, loc='upper right')

    plt.tight_layout()
    plt.savefig('phonon_04_quantum_field.png', dpi=300)
    print("Saved figure: phonon_04_quantum_field.png")

if __name__ == '__main__':
    plot_quantum_fluctuations()
