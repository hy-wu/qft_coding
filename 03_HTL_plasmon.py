"""
HTL Self-Energy and Plasmon Dispersion in QCD
==============================================
From scalar tadpole to gauge theory:
the Hard Thermal Loop (HTL) gluon self-energy.

Physics:
  At high T, the QCD medium modifies gluon propagation:
  - Longitudinal (electric) mode: Debye screening, m_D ~ gT
  - Transverse (magnetic) mode: dynamical screening, plasma frequency

  Both modes satisfy different dispersion relations and merge
  at zero momentum into the plasma frequency ω_pl = m_D/√3.

References:
  [1] Kapusta & Gale Sec 5.5
  [2] Le Bellac Sec 4.3
  [3] Laine & Vuorinen Sec 4.3
"""

import numpy as np
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

# ================================================================
# Parameters: QCD with N_c=3, N_f=3
# ================================================================
g = 2.0             # coupling, roughly alpha_s ~ 0.3
T = 1.0             # temperature (scale setter)
Nc = 3              # colors
Nf = 3              # active flavors

# Debye mass squared
# m_D^2 = g^2 T^2 (N_c/3 + N_f/6)
mD2 = g**2 * T**2 * (Nc / 3.0 + Nf / 6.0)
mD = np.sqrt(mD2)

print(f"m_D = {mD:.4f}")
print(f"w_pl = m_D/rt3 = {mD/np.sqrt(3):.4f}")

# ================================================================
# HTL self-energies
#
# For |ω| > p (propagating modes): real analytic
# For |ω| < p (Landau damped):   develop imaginary part
#
# We only need the real part for propagating mode dispersion.
# ================================================================

def Pi_L(omega, p):
    """
    Longitudinal (electric) HTL self-energy.
    Π_L(ω,p) = m_D^2 [1 - ω/(2p) ln((ω+p)/(ω-p))]

    For ω > p: real
    For ω < p: complex (→ Landau damping)
    """
    if p < 1e-15:
        # p -> 0 limit: Π_L(ω,0) = m_D^2 * (ω^2)/(3*ω^2) ... actually
        # the p=0 limit is subtle; we handle it via the small-p expansion
        return mD2 / 3.0

    if omega <= p:
        # In the light-cone or spacelike region → complex
        # The imaginary part signals Landau damping
        # We only use the real part for propagating mode root-finding
        return np.nan

    # ω > p: real propagating mode
    log_val = np.log((omega + p) / (omega - p))
    return mD2 * (1.0 - omega / (2.0 * p) * log_val)


def Pi_T(omega, p):
    """
    Transverse (magnetic) HTL self-energy.

    Π_T(ω,p) = m_D^2/2 [ ω²/p² - ω(ω²-p²)/(2p³) ln((ω+p)/(ω-p)) ]

    For ω > p: real
    For ω < p: complex
    """
    if p < 1e-15:
        return mD2 / 3.0

    if omega <= p:
        return np.nan

    log_val = np.log((omega + p) / (omega - p))

    term1 = omega**2 / p**2
    term2 = omega * (omega**2 - p**2) / (2.0 * p**3) * log_val

    return 0.5 * mD2 * (term1 - term2)


# ================================================================
# Dispersion relations
#
# Longitudinal (plasmon):
#   Δ_L^{-1} = p^2 + Π_L(ω, p) = 0   →   ω = ω_L(p)
#
# Transverse:
#   Δ_T^{-1} = -ω^2 + p^2 + Π_T(ω, p) = 0   →   ω = ω_T(p)
# ================================================================

def f_L(omega, p):
    """Longitudinal dispersion function: zero = plasmon mode"""
    PiL = Pi_L(omega, p)
    if np.isnan(PiL):
        return np.nan
    return p**2 + PiL


def f_T(omega, p):
    """Transverse dispersion function: zero = transverse mode"""
    PiT = Pi_T(omega, p)
    if np.isnan(PiT):
        return np.nan
    return -omega**2 + p**2 + PiT


# ================================================================
# Root finding
#
# Longitudinal mode exists for ω > p.
# At ω = p⁺: Π_L → -∞   →   f_L → -∞
# At ω → ∞:  Π_L → 0    →   f_L → p² > 0
# So there is exactly one root in (p, ∞).
#
# Transverse mode:
# At ω = p⁺: Π_T → m_D²/2   →   f_T → m_D²/2 > 0
# At ω → ∞:  Π_T → 0        →   f_T → -ω² → -∞
# So there is exactly one root in (p, ∞).
# ================================================================

def solve_dispersion(func, p, omega_guess=None):
    """
    Solve dispersion relation for a given momentum p.
    Uses bisection bracketing.
    """
    # Find lower bound: function must be negative here (for longitudinal)
    # or positive (for transverse), then go to where it flips sign.

    # Start very close to p from above
    lo = p * 1.0001 + 1e-6
    f_lo = func(lo, p)

    # Find upper bound where function has opposite sign
    hi = lo + 0.5
    for _ in range(30):
        f_hi = func(hi, p)
        if np.isnan(f_hi):
            hi *= 1.5
            continue
        if f_lo * f_hi < 0:
            break
        hi *= 1.5

    if f_lo * f_hi > 0:
        # Failed to bracket — maybe this p has no propagating mode
        return np.nan

    try:
        sol = root_scalar(func, args=(p,), bracket=(lo, hi),
                          method='bisect', xtol=1e-10)
        return sol.root
    except (ValueError, RuntimeError):
        return np.nan


# ================================================================
# Scan over momenta
# ================================================================
momentum_grid = np.linspace(0.0, 5.0 * mD, 200)

omega_L = []
omega_T = []

# p=0 is special: both modes give ω = m_D/√3
omega_0 = mD / np.sqrt(3)
omega_L.append(omega_0)
omega_T.append(omega_0)

# p>0: solve dispersion relations
for p in momentum_grid[1:]:
    # Transverse is easier: f_T(p⁺) > 0, f_T(∞) < 0
    wT = solve_dispersion(f_T, p)
    omega_T.append(wT)

    # Longitudinal: f_L(p⁺) < 0, f_L(∞) > 0
    wL = solve_dispersion(f_L, p)
    omega_L.append(wL)

    if np.isfinite(wL) and np.isfinite(wT):
        print(f"p={p:.3f}:  wL={wL:.4f}  wT={wT:.4f}")

omega_L = np.array(omega_L)
omega_T = np.array(omega_T)

# ================================================================
# Asymptotic forms
# ================================================================

# Large-p asymptotic for transverse mode:
# ω_T² ≈ p² + m_D²/3  (simple screening approximation)
omega_T_asym = np.sqrt(momentum_grid**2 + mD2 / 3.0)

# Light cone
light_cone = momentum_grid.copy()
light_cone[0] = 1e-15  # avoid 0/0 in ratio plots

# ================================================================
# Plots
# ================================================================

# --- Plot 1: Dispersion relations ---
plt.figure(figsize=(8, 6))

plt.plot(momentum_grid, omega_L, '-', color='#1f77b4', linewidth=2.5,
         label=r'Longitudinal (plasmon)')
plt.plot(momentum_grid, omega_T, '-', color='#d62728', linewidth=2.5,
         label=r'Transverse')
plt.plot(momentum_grid, omega_T_asym, '--', color='#2ca02c', linewidth=1.5,
         label=r'$\sqrt{p^2 + m_D^2/3}$')
plt.plot(momentum_grid, light_cone, ':', color='gray', linewidth=1.5,
         label=r'Light cone $\omega = p$')

# Mark the plasma frequency
plt.axhline(y=omega_0, color='purple', linestyle='--', alpha=0.5)
plt.annotate(r'$\omega_{\mathrm{pl}} = m_D/\sqrt{3}$',
             xy=(0.1, omega_0 + 0.05*mD), fontsize=12, color='purple')

plt.xlabel(r'$p$', fontsize=14)
plt.ylabel(r'$\omega(p)$', fontsize=14)
plt.title(r'HTL Plasmon Dispersion in QCD ($N_c=3, N_f=3$)', fontsize=13)
plt.legend(fontsize=11)
plt.xlim(0, momentum_grid[-1])
plt.ylim(0, momentum_grid[-1] * 1.1)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('03_HTL_dispersion.png', dpi=150)

# --- Plot 2: Approach to light cone ---
plt.figure(figsize=(8, 5))

# Only plot regions where we have valid solutions
mask_L = np.isfinite(omega_L)
mask_T = np.isfinite(omega_T)

delta_L = omega_L - light_cone
delta_T = omega_T - light_cone

plt.semilogy(momentum_grid[mask_L], delta_L[mask_L], '-',
             color='#1f77b4', linewidth=2, label=r'$\omega_L - p$')
plt.semilogy(momentum_grid[mask_T], delta_T[mask_T], '-',
             color='#d62728', linewidth=2, label=r'$\omega_T - p$')

plt.axvline(x=mD, color='gray', linestyle=':', alpha=0.5, label=r'$m_D$')
plt.xlabel(r'$p$', fontsize=14)
plt.ylabel(r'$\omega(p) - p$', fontsize=14)
plt.title(r'Deviation from Light Cone', fontsize=13)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('03_HTL_asymptotic.png', dpi=150)

# --- Plot 3: Thermal mass comparison ---
# Compare scalar thermal mass from Ex01 with gluonic Debye mass
plt.figure(figsize=(7, 4))

# The scalar thermal mass from φ⁴: m_th² = λT²/24
# For comparison, scale so they both appear on same plot
scalar_mth2 = 1.0 * T**2 / 24.0
qcd_mD2_vals = [g**2 * T**2 * (3/3 + f/6) for f in [0, 2, 3, 6]]
labels = [r'$N_f=0$ (pure glue)', r'$N_f=2$', r'$N_f=3$', r'$N_f=6$']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

x = np.arange(len(qcd_mD2_vals))
width = 0.35

plt.bar(x, [v/T**2 for v in qcd_mD2_vals], width, color=colors, alpha=0.8)
plt.axhline(y=scalar_mth2/T**2, color='purple', linestyle='--', linewidth=2,
            label=r'$\phi^4$: $\lambda T^2/24$')

plt.xticks(x, labels, fontsize=11)
plt.ylabel(r'$m_D^2 / T^2$', fontsize=14)
plt.title(r'Debye Mass in QCD vs Scalar Thermal Mass', fontsize=13)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig('03_mD_comparison.png', dpi=150)

print()
print("Saved: 03_HTL_dispersion.png")
print("Saved: 03_HTL_asymptotic.png")
print("Saved: 03_mD_comparison.png")
print()
print("--- Summary ---")
print(f"Debye mass:         m_D = {mD:.3f}")
print(f"Plasma frequency:   w_pl = m_D/rt3 = {omega_0:.3f}")
print(f"At p=0: both modes degenerate at w_pl")
print(f"Asymptotically (large p):")
print(f"  Longitudinal:  w_L ~ p (exponentially screened)")
print(f"  Transverse:    w_T ~ p (w_T^2 ~ p^2 + m_D^2/3)")
