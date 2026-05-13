"""
Scalar Spectral Function at Finite Temperature
===============================================
Extend the tadpole result to compute the full spectral function.

Physics flow:
  1. Tadpole gives thermal mass m_th^2(T) = λT²/24 (high-T)
  2. Beyond tadpole: self-energy gets an imaginary part from
     2->2 scattering at finite T → thermal width Γ(T)
  3. Full propagator: D_R⁻¹ = -ω² + p² + m₀² + Π_R(ω,p)
  4. Spectral function: ρ(ω,p) = -(1/π) Im D_R(ω,p)
  5. As T increases: peak shifts (mass) AND broadens (width)

Connection to QGP:
  Exactly analogous to quarkonium (J/ψ, Υ) in QGP:
  - Colored screening → peak shift (mass modification)
  - Gluon dissociation → peak broadening (melting)
  Experimental signature: dilepton spectrum shows broadened ρ peak
"""

import numpy as np
import matplotlib.pyplot as plt

# ================================================================
# Parameters (natural units)
# ================================================================
LAMBDA = 1.0          # coupling
M0 = 0.5              # bare mass

# ================================================================
# Thermal mass from tadpole (Exercise 1)
# High-T limit: m_th^2 = λT²/24
# ================================================================
def thermal_mass_sq(T):
    """Thermal mass squared from tadpole self-energy"""
    if T <= 1e-10:
        return 0.0
    return LAMBDA * T**2 / 24.0

# ================================================================
# Thermal width:
# In φ⁴ theory at high T, the damping rate (width) is
#   Γ ∼ λ² T / (4π)    (leading-order from 2->2 scattering)
#
# The full self-energy near the quasi-particle pole:
#   Π_R(ω, p) ≈ m_th² − i ω Γ
#
# This gives the Breit-Wigner spectral function:
#   ρ(ω,p) = 1/π · ωΓ / [(ω²−ω_p²)² + (ωΓ)²]
#   where ω_p² = p² + M0² + m_th²(T)
# ================================================================
def thermal_width(T, g=LAMBDA):
    """Leading-order thermal width for scalar φ⁴"""
    if T <= 1e-10:
        return 0.0
    return g**2 * T / (4.0 * np.pi)


# ================================================================
# Spectral function ρ(ω, p; T)
#
# Retarded propagator:
#   D_R(ω, p) = 1 / (−ω² + p² + M0² + m_th² − i ω Γ)
#
# Spectral function:
#   ρ(ω, p) = −(1/π) Im D_R(ω, p)
#            = 1/π · ωΓ / [(ω² − ω_p²)² + (ωΓ)²]
#
# Normalization: ∫_{-∞}^{∞} dω ρ(ω, p) = Z(p)  (波函數重整化)
# In free limit Γ→0+: ρ(ω,p) = sgn(ω) δ(ω²−ω_p²)
# ================================================================
def spectral_function(omega, p, T):
    """
    Scalar spectral function at finite temperature.

    Parameters:
        omega: real frequency
        p: spatial momentum magnitude
        T: temperature
    """
    mth2 = thermal_mass_sq(T)
    Gamma = thermal_width(T)

    # quasi-particle dispersion
    omega_p_sq = p**2 + M0**2 + mth2

    # Breit-Wigner form
    denom = (omega**2 - omega_p_sq)**2 + (omega * Gamma)**2
    rho = (omega * Gamma) / (np.pi * denom)

    return rho


# ================================================================
# Scan: fixed p, vary T
# ================================================================
p_fixed = 0.0                       # zero momentum (rest frame)
omega_grid = np.linspace(0, 4.0, 400)

T_values = [0.1, 0.5, 1.0, 2.0, 4.0]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

plt.figure(figsize=(8, 5))

for T, color in zip(T_values, colors):
    rho = np.array([spectral_function(w, p_fixed, T) for w in omega_grid])
    plt.plot(omega_grid, rho, '-', color=color, linewidth=2,
             label=rf'$T={T:.1f},\; m_{{\rm th}}={np.sqrt(thermal_mass_sq(T)):.2f}$')

# Free theory line (T=0 limit: delta function at ω = ω₀)
# We approximate by a very narrow Breit-Wigner
T0 = 1e-3
rho_free = np.array([spectral_function(w, p_fixed, T0) for w in omega_grid])
# scale free part for visibility
# rho_free *= 5

plt.axvline(x=np.sqrt(M0**2), color='gray', linestyle=':', alpha=0.5,
            label=r'$\omega = m_0$ (free)')

plt.xlabel(r'$\omega$', fontsize=14)
plt.ylabel(r'$\rho(\omega, p=0)$', fontsize=14)
plt.title(r'Scalar Spectral Function $T$ dependence', fontsize=13)
plt.legend(fontsize=11)
plt.ylim(0, None)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('02_spectral_function_T.png', dpi=150)


# ================================================================
# Scan: fixed T, vary p  (plasmon dispersion)
# ================================================================
T_fixed = 2.0
momenta = [0.0, 1.0, 2.0, 4.0]

plt.figure(figsize=(8, 5))

for p, color in zip(momenta, colors):
    rho = np.array([spectral_function(w, p, T_fixed) for w in omega_grid])
    omega_p = np.sqrt(p**2 + M0**2 + thermal_mass_sq(T_fixed))
    plt.plot(omega_grid, rho, '-', color=color, linewidth=2,
             label=rf'$p={p:.1f},\; \omega_p={omega_p:.2f}$')

plt.axvline(x=0, color='gray', linestyle=':', alpha=0.3)
plt.xlabel(r'$\omega$', fontsize=14)
plt.ylabel(r'$\rho(\omega, p)$', fontsize=14)
plt.title(rf'Spectral Function at $T={T_fixed}$', fontsize=13)
plt.legend(fontsize=11)
plt.ylim(0, None)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('02_spectral_function_p.png', dpi=150)


# ================================================================
# Key physics: peak position vs width as function of T
# ================================================================
T_range = np.linspace(0.1, 5.0, 50)
peak_positions = []
peak_widths = []

for T in T_range:
    omega_p = np.sqrt(M0**2 + thermal_mass_sq(T))
    Gamma = thermal_width(T)
    peak_positions.append(omega_p)
    peak_widths.append(Gamma)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(T_range, peak_positions, 'o-', color='#1f77b4', linewidth=2, markersize=4)
ax1.set_xlabel(r'$T$', fontsize=14)
ax1.set_ylabel(r'$\omega_p$ (peak position)', fontsize=14)
ax1.set_title(r'Thermal mass shift', fontsize=13)
ax1.grid(True, alpha=0.3)

ax2.plot(T_range, peak_widths, 's-', color='#d62728', linewidth=2, markersize=4)
ax2.set_xlabel(r'$T$', fontsize=14)
ax2.set_ylabel(r'$\Gamma$ (FWHM)', fontsize=14)
ax2.set_title(r'Thermal broadening $\sim\lambda^2 T$', fontsize=13)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('02_peak_T_dependence.png', dpi=150)

print("Saved: 02_spectral_function_T.png")
print("Saved: 02_spectral_function_p.png")
print("Saved: 02_peak_T_dependence.png")
print()
print("--- Summary ---")
print("Free particle: delta-function peak at omega = sqrt(p^2 + m0^2)")
print("At T>0: peak shifts to omega = sqrt(p^2 + m0^2 + lambda T^2/24)")
print("        peak broadens with width Gamma = lambda^2 T/(4pi)")
print()
print("This is the scalar analogue of quarkonium melting in QGP:")
print("  J/psi, Upsilon peaks broaden and disappear above T_c")
print("  Measured in dilepton channel at LHC/RHIC")
