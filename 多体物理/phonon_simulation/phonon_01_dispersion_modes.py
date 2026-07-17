r"""
Phonon Theory Script 01: 1D Monoatomic and Diatomic Lattice Dispersion & Normal Modes
===================================================================================
Formula to Code Mapping:
- Eq (1.1) 1D Monoatomic Dispersion: \omega(k) = 2 * sqrt(C/M) * |sin(k*a/2)|
  -> Function `dispersion_monoatomic(k, C, M, a)`
- Eq (1.2) 1D Diatomic Dynamical Matrix:
  D(k) = [[ 2C/M1, -C/M1*(1+exp(-i*k*a)) ],
          [ -C/M2*(1+exp(i*k*a)), 2C/M2 ]]
  -> Function `dynamical_matrix_diatomic(k, C, M1, M2, a)`
- Eq (1.3) Normal Modes Displacement:
  u_n(t) = Re[ e_{k,\sigma} * exp(i*(k*R_n - \omega*t)) ]
  -> Function `simulate_atomic_motion(...)`
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150

def dispersion_monoatomic(k, C=1.0, M=1.0, a=1.0):
    return 2.0 * np.sqrt(C / M) * np.abs(np.sin(k * a / 2.0))

def dispersion_diatomic(k, C=1.0, M1=1.0, M2=2.5, a=1.0):
    k = np.atleast_1d(k)
    omega_opt = np.zeros_like(k, dtype=float)
    omega_ac  = np.zeros_like(k, dtype=float)
    eigvecs_opt = np.zeros((len(k), 2), dtype=complex)
    eigvecs_ac  = np.zeros((len(k), 2), dtype=complex)

    for i, ki in enumerate(k):
        D = np.array([
            [2.0 * C / M1, -C / M1 * (1.0 + np.exp(-1j * ki * a))],
            [-C / M2 * (1.0 + np.exp(1j * ki * a)), 2.0 * C / M2]
        ], dtype=complex)

        eigvals, eigvecs = np.linalg.eigh(D)
        idx = np.argsort(eigvals)
        omega_ac[i] = np.sqrt(np.maximum(0.0, eigvals[idx[0]]))
        omega_opt[i] = np.sqrt(np.maximum(0.0, eigvals[idx[1]]))
        
        eigvecs_ac[i] = eigvecs[:, idx[0]]
        eigvecs_opt[i] = eigvecs[:, idx[1]]

    return omega_opt, omega_ac, eigvecs_opt, eigvecs_ac

def plot_dispersion_and_modes():
    a = 1.0
    C = 1.0
    M1, M2 = 1.0, 2.5
    k_vals = np.linspace(-np.pi / a, np.pi / a, 500)

    omega_mono = dispersion_monoatomic(k_vals, C=C, M=M1, a=a)
    omega_opt, omega_ac, vecs_opt, vecs_ac = dispersion_diatomic(k_vals, C=C, M1=M1, M2=M2, a=a)

    fig = plt.figure(figsize=(14, 6))

    # --- Panel 1: Dispersion Curves ---
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(k_vals / (np.pi / a), omega_mono, 'k--', label=r'单原子链 $\omega(k)$ ($M=1.0$)', alpha=0.7, lw=2)
    ax1.plot(k_vals / (np.pi / a), omega_ac, 'b-', label=r'双原子链 声学支 $\omega_-(k)$', lw=2.5)
    ax1.plot(k_vals / (np.pi / a), omega_opt, 'r-', label=r'双原子链 光学支 $\omega_+(k)$', lw=2.5)

    w_ac_max = np.sqrt(2 * C / M2)
    w_opt_min = np.sqrt(2 * C / M1)
    w_opt_max = np.sqrt(2 * C * (1/M1 + 1/M2))

    ax1.axhline(w_ac_max, color='blue', linestyle=':', alpha=0.6)
    ax1.axhline(w_opt_min, color='red', linestyle=':', alpha=0.6)
    ax1.fill_between(k_vals / (np.pi / a), w_ac_max, w_opt_min, color='gray', alpha=0.15, label=r'声子能隙 (Band Gap)')

    ax1.set_xlabel(r'波矢 $k / (\pi / a)$ [第一布里渊区]', fontsize=12)
    ax1.set_ylabel(r'角频率 $\omega$', fontsize=12)
    ax1.set_title(r'图 1a: 1D 晶格声子色散曲线 $\omega(k)$', fontsize=14, pad=12)
    ax1.set_xlim([-1, 1])
    ax1.set_ylim([0, w_opt_max * 1.08])
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend(fontsize=10, loc='upper center')

    vs = a * np.sqrt(C / (M1 + M2))
    ax1.annotate(f'长波极限声速 $v_s = \\frac{{d\\omega}}{{dk}} = {vs:.3f}$',
                 xy=(0.0, 0.0), xytext=(0.1, 0.3),
                 arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=6),
                 fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="yellow", ec="b", lw=1, alpha=0.5))

    # --- Panel 2: Real-space Atomic Motion ---
    ax2 = fig.add_subplot(1, 2, 2)
    N_cells = 6
    x_M1 = np.arange(N_cells) * a
    x_M2 = x_M1 + 0.5 * a

    k_demo = 0.25 * np.pi / a
    _, _, vecs_opt_demo, vecs_ac_demo = dispersion_diatomic(k_demo, C=C, M1=M1, M2=M2, a=a)
    
    e_ac = vecs_ac_demo[0]
    e_opt = vecs_opt_demo[0]

    u_ac_M1 = np.real(e_ac[0] * np.exp(1j * k_demo * x_M1))
    u_ac_M2 = np.real(e_ac[1] * np.exp(1j * k_demo * x_M2))

    u_opt_M1 = np.real(e_opt[0] * np.exp(1j * k_demo * x_M1))
    u_opt_M2 = np.real(e_opt[1] * np.exp(1j * k_demo * x_M2))

    u_ac_M1 /= np.max(np.abs(u_ac_M1)) * 2.5
    u_ac_M2 /= np.max(np.abs(u_ac_M2)) * 2.5
    u_opt_M1 /= np.max(np.abs(u_opt_M1)) * 2.5
    u_opt_M2 /= np.max(np.abs(u_opt_M2)) * 2.5

    y_ac = 1.0
    y_opt = 0.0

    ax2.scatter(x_M1, np.full_like(x_M1, y_ac), s=200, color='royalblue', label=r'原子 $M_1$', zorder=3)
    ax2.scatter(x_M2, np.full_like(x_M2, y_ac), s=350, color='crimson', label=r'原子 $M_2$', zorder=3)
    
    # Quiver fixed with matching y vector
    ax2.quiver(x_M1, np.full_like(x_M1, y_ac), u_ac_M1, np.zeros_like(u_ac_M1), angles='xy', scale_units='xy', scale=1, color='blue', width=0.008)
    ax2.quiver(x_M2, np.full_like(x_M2, y_ac), u_ac_M2, np.zeros_like(u_ac_M2), angles='xy', scale_units='xy', scale=1, color='red', width=0.008)

    ax2.scatter(x_M1, np.full_like(x_M1, y_opt), s=200, color='royalblue', zorder=3)
    ax2.scatter(x_M2, np.full_like(x_M2, y_opt), s=350, color='crimson', zorder=3)
    ax2.quiver(x_M1, np.full_like(x_M1, y_opt), u_opt_M1, np.zeros_like(u_opt_M1), angles='xy', scale_units='xy', scale=1, color='blue', width=0.008)
    ax2.quiver(x_M2, np.full_like(x_M2, y_opt), u_opt_M2, np.zeros_like(u_opt_M2), angles='xy', scale_units='xy', scale=1, color='red', width=0.008)

    ax2.text(-0.5, y_ac, "声学模 (同向运动):\n$M_1, M_2$ 同相位", fontsize=11, fontweight='bold', va='center', ha='right')
    ax2.text(-0.5, y_opt, "光学模 (反向运动):\n$M_1, M_2$ 反相位", fontsize=11, fontweight='bold', va='center', ha='right')

    ax2.set_title(r'图 1b: 极化矢量与实空间振动图像 ($k \to 0$)', fontsize=14, pad=12)
    ax2.set_xlabel(r'晶格位置 $x / a$', fontsize=12)
    ax2.set_yticks([])
    ax2.set_xlim([-2.0, N_cells * a])
    ax2.set_ylim([-0.8, 1.8])
    ax2.grid(True, linestyle=':', alpha=0.4)
    ax2.legend(fontsize=10, loc='upper right')

    plt.tight_layout()
    plt.savefig('phonon_01_dispersion_modes.png', dpi=300)
    print("Saved figure: phonon_01_dispersion_modes.png")

if __name__ == '__main__':
    plot_dispersion_and_modes()
