"""
Phonon Theory Script 02: Phonon Density of States (DOS) and Van Hove Singularities
================================================================================
Formula to Code Mapping:
- Eq (2.1) General DOS definition:
  g(\omega) = \int \frac{d^d k}{(2\pi)^d} \delta(\omega - \omega(k)) = \int_{\omega(k)=\omega} \frac{dS}{(2\pi)^d |\nabla_k \omega(k)|}
- Eq (2.2) 1D Analytical DOS:
  g_{1D}(\omega) = \frac{2}{\pi} \frac{1}{\sqrt{\omega_{max}^2 - \omega^2}}
  -> Function `dos_1d_analytical(\omega, \omega_max)`
- Eq (2.3) 2D Numerical DOS Sampling on Brillouin Zone (Square Lattice):
  \omega(k_x, k_y) = \omega_0 \sqrt{\sin^2(k_x a/2) + \sin^2(k_y a/2)}
  -> Function `dos_2d_numerical(...)`
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150

def dos_1d_analytical(omega, omega_max=2.0):
    """
    Formula (2.2): Analytical 1D DOS with 1/\sqrt{\omega_{max}^2 - \omega^2} singularity.
    """
    g = np.zeros_like(omega, dtype=float)
    mask = (omega >= 0) & (omega < omega_max)
    g[mask] = (2.0 / np.pi) / np.sqrt(omega_max**2 - omega[mask]**2)
    return g

def compute_2d_dispersion_and_dos(N_k=1000, n_bins=300):
    """
    Formula (2.3): 2D Square Lattice Phonon Dispersion:
    \omega(k_x, k_y) = \omega_{max} \sqrt{ (\sin^2(k_x a/2) + \sin^2(k_y a/2)) / 2 }
    """
    omega_max = 2.0
    kx = np.linspace(0, np.pi, N_k)
    ky = np.linspace(0, np.pi, N_k)
    Kx, Ky = np.meshgrid(kx, ky)

    # 2D dispersion relation
    Omega_2D = omega_max * np.sqrt(0.5 * (np.sin(Kx / 2.0)**2 + np.sin(Ky / 2.0)**2))

    # Numerical Histogram for DOS g(\omega)
    hist, bin_edges = np.histogram(Omega_2D, bins=n_bins, range=(0, omega_max), density=True)
    omega_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    return Kx, Ky, Omega_2D, omega_centers, hist

def plot_dos_and_van_hove():
    fig = plt.figure(figsize=(15, 5))

    # --- Panel 1: 1D Dispersion & 1D DOS (1/sqrt singularity) ---
    ax1 = fig.add_subplot(1, 3, 1)
    w_max = 2.0
    k_1d = np.linspace(0, np.pi, 300)
    w_1d = w_max * np.sin(k_1d / 2.0)
    
    omega_grid = np.linspace(0, w_max * 0.999, 500)
    g_1d = dos_1d_analytical(omega_grid, w_max)

    ax1.plot(g_1d, omega_grid, 'r-', lw=2.5, label='1D 声子态密度 $g(\omega)$')
    ax1.axhline(w_max, color='black', linestyle='--', label=r'$\omega_{max} (\nabla_k \omega = 0)$')
    ax1.set_xlabel('态密度 $g(\omega)$', fontsize=12)
    ax1.set_ylabel('角频率 $\omega$', fontsize=12)
    ax1.set_title('图 2a: 1D 态密度与 Van Hove 奇点\n($g(\omega) \sim (\omega_{max}-\omega)^{-1/2}$)', fontsize=13)
    ax1.set_xlim([0, 3.5])
    ax1.set_ylim([0, w_max * 1.1])
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend(fontsize=10, loc='lower right')

    # --- Panel 2: 2D Dispersion Contour Surface ---
    Kx, Ky, Omega_2D, omega_centers, g_2d = compute_2d_dispersion_and_dos(N_k=800, n_bins=250)
    
    ax2 = fig.add_subplot(1, 3, 2)
    contour = ax2.contourf(Kx / np.pi, Ky / np.pi, Omega_2D, levels=25, cmap='viridis')
    cbar = plt.colorbar(contour, ax=ax2)
    cbar.set_label('频率 $\omega$', fontsize=11)
    
    # Mark critical points (Van Hove saddle point and extremum)
    ax2.scatter([1.0], [0.0], color='red', s=100, zorder=5, label='鞍点 $(\pi, 0)$ [Van Hove奇点]')
    ax2.scatter([0.0], [1.0], color='red', s=100, zorder=5)
    ax2.scatter([1.0], [1.0], color='cyan', s=100, zorder=5, label='极大点 $(\pi, \pi)$')
    
    ax2.set_xlabel('$k_x / \pi$', fontsize=12)
    ax2.set_ylabel('$k_y / \pi$', fontsize=12)
    ax2.set_title('图 2b: 2D 二维网格色散等频线 $\omega(k_x, k_y)$', fontsize=13)
    ax2.legend(fontsize=9, loc='upper left')

    # --- Panel 3: 2D DOS with Logarithmic Singularity ---
    ax3 = fig.add_subplot(1, 3, 3)
    ax3.plot(g_2d, omega_centers, 'b-', lw=2.5, label='2D 数值态密度 $g(\omega)$')
    
    # Saddle point frequency is \omega_saddle = \omega_max * sqrt(0.5 * (1 + 0)) = \omega_max / sqrt(2) = 2.0 / 1.414 = 1.414
    w_saddle = w_max / np.sqrt(2.0)
    ax3.axhline(w_saddle, color='crimson', linestyle='--', label=r'鞍点频率 $\omega_{saddle} = \sqrt{2}$')
    ax3.axhline(w_max, color='black', linestyle=':', label=r'带顶频率 $\omega_{max} = 2.0$')

    ax3.annotate('对数奇点\n logarithmic singularity\n($\\nabla_k \\omega = 0$)',
                 xy=(g_2d[np.argmin(np.abs(omega_centers - w_saddle))], w_saddle),
                 xytext=(1.2, 1.1),
                 arrowprops=dict(facecolor='red', shrink=0.08, width=1.5, headwidth=6),
                 fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="pink", ec="r", lw=1, alpha=0.6))

    ax3.set_xlabel('态密度 $g(\omega)$', fontsize=12)
    ax3.set_ylabel('角频率 $\omega$', fontsize=12)
    ax3.set_title('图 2c: 2D 态密度中 Van Hove 奇点', fontsize=13)
    ax3.set_xlim([0, np.max(g_2d) * 1.15])
    ax3.set_ylim([0, w_max * 1.1])
    ax3.grid(True, linestyle='--', alpha=0.5)
    ax3.legend(fontsize=10, loc='lower right')

    plt.tight_layout()
    plt.savefig('phonon_02_dos_van_hove.png', dpi=300)
    print("Saved figure: phonon_02_dos_van_hove.png")

if __name__ == '__main__':
    plot_dos_and_van_hove()
