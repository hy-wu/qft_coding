"""
Phase 2.2: Lattice phi^4 — The Equal-Time Propagator (2-point correlator)
===========================================================================

Physics:
  C_{ij} = <0| phi_i phi_j |0>    (equal-time ground state correlator)

  Free theory:
    C^0_{ij} = sum_k U_{ik} U_{jk} / (2 omega_k)
             ~ exp(-m * |i-j|) at large separation
    Correlation length: xi = 1/m

  Interacting theory (lambda > 0):
    C_{ij} = <psi_0| phi_i phi_j |psi_0>
    Mass renormalization: the physical mass changes due to lambda phi^4
    This is the propagator — the building block of all Feynman diagrams!

Connection to QFT:
  In the continuum:  D(x-y) = int d^3k/(2pi)^3 e^{ik.(x-y)} / (2 omega_k)
  On the lattice:    C_{ij} = (1/N) sum_k cos(k.(i-j)) / (2 omega_k)
"""

import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# ================================================================
# Parameters
# ================================================================
m = 0.5                     # bare mass
g = 1.0                     # nearest-neighbor coupling
lam = 0.2                   # phi^4 coupling (larger than 06 to see effect)
N_sites = 8                 # lattice sites
N_fock = 3                  # Fock truncation per mode (3^8 = 6561)

# ================================================================
# Step 1: Free normal modes
# ================================================================
K = np.zeros((N_sites, N_sites))
for i in range(N_sites):
    K[i, i] = m**2 + 2 * g
    K[i, (i+1) % N_sites] = -g
    K[i, (i-1) % N_sites] = -g

omega2_k, U = eigh(K)
omega_k = np.sqrt(omega2_k)
idx_sort = np.argsort(omega_k)
omega_k = omega_k[idx_sort]
U = U[:, idx_sort]

print("=== Phase 2.2: Lattice phi^4 propagator ===")
print(f"N_sites = {N_sites},  m = {m},  g = {g},  lambda = {lam}")
print(f"Mode frequencies:")
for k in range(N_sites):
    print(f"  omega_{k} = {omega_k[k]:.4f}")

# ================================================================
# Step 2: Free equal-time correlator
# C^0_{ij} = sum_k U_{ik} U_{jk} / (2 omega_k)
#
# This IS the free propagator at equal time:
#   D(x_i - x_j, t=0) = <0| phi_i phi_j |0>
# ================================================================
C_free = np.zeros((N_sites, N_sites))
for k in range(N_sites):
    C_free += np.outer(U[:, k], U[:, k]) / (2.0 * omega_k[k])

# Check that it's symmetric and positive-definite
print(f"\nFree correlator C^0_{0}{0} = {C_free[0,0]:.6f}")
print(f"C^0 symmetric: {np.allclose(C_free, C_free.T)}")

# ================================================================
# Step 3: Build interacting Hamiltonian
# ================================================================
dim = N_fock ** N_sites
print(f"\nHilbert space dimension: {dim}")

def idx_to_state(idx):
    state = []
    for _ in range(N_sites):
        state.append(idx % N_fock)
        idx //= N_fock
    return tuple(state)

# H_0 diagonal
H0_diag = np.zeros(dim)
for idx in range(dim):
    state = idx_to_state(idx)
    E = sum(omega_k[k] * (n + 0.5) for k, n in enumerate(state))
    H0_diag[idx] = E

# Single-mode x matrix
def build_x_matrix(Nf, omega):
    x = np.zeros((Nf, Nf))
    for n in range(Nf):
        if n > 0:
            x[n, n-1] = np.sqrt(n)
        if n < Nf - 1:
            x[n, n+1] = np.sqrt(n + 1)
    return x / np.sqrt(2 * omega)

# X_k operators in full space
X_ops = []
for k in range(N_sites):
    x_k_single = build_x_matrix(N_fock, omega_k[k])
    op_list = [np.eye(N_fock)] * N_sites
    op_list[k] = x_k_single
    X_k = op_list[0]
    for op in op_list[1:]:
        X_k = np.kron(X_k, op)
    X_ops.append(X_k)

# phi_i operators
phi_ops = []
for i in range(N_sites):
    phi_i = sum(U[i, k] * X_ops[k] for k in range(N_sites))
    phi_ops.append(phi_i)

# H_1 = lambda/4! sum_i phi_i^4
H1 = np.zeros((dim, dim))
for i in range(N_sites):
    phi2 = phi_ops[i] @ phi_ops[i]
    phi4 = phi2 @ phi2
    H1 += phi4
H1 *= lam / 24.0

# Full Hamiltonian
H_full = np.diag(H0_diag) + H1

print("Diagonalizing...")
E_exact, psi_exact = eigh(H_full)
psi0 = psi_exact[:, 0]
E0 = E_exact[0]
print(f"Exact GS energy: E0 = {E0:.6f}")

# ================================================================
# Step 4: Interacting correlator
# C_{ij} = <psi_0| phi_i phi_j |psi_0>
# ================================================================
C_exact = np.zeros((N_sites, N_sites))
for i in range(N_sites):
    phi_i_psi = phi_ops[i] @ psi0
    for j in range(N_sites):
        C_exact[i, j] = np.dot(phi_i_psi.conj(), phi_ops[j] @ psi0)

# Also compute <phi_i^2> for each site (check translation invariance)
phi2_exact = np.array([C_exact[i, i] for i in range(N_sites)])
phi2_free = np.array([C_free[i, i] for i in range(N_sites)])

print(f"\n<phi_i^2> (free):     mean = {phi2_free.mean():.6f}, "
      f"spread = {np.ptp(phi2_free):.2e}")
print(f"<phi_i^2> (exact):    mean = {phi2_exact.mean():.6f}, "
      f"spread = {np.ptp(phi2_exact):.2e}")

# Analytic free result: 1/N sum_k 1/(2 omega_k)
phi2_free_analytic = sum(1.0/(2*w) for w in omega_k) / N_sites
print(f"<phi_i^2> (analytic): {phi2_free_analytic:.6f}")

# ================================================================
# Step 5: Plot the correlator C(r) = C_{0, r}
# ================================================================
r_values = np.arange(N_sites)
C_free_r = C_free[0, :]
C_exact_r = C_exact[0, :]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Left: linear scale
ax1.plot(r_values, C_free_r, 'o-', color='#1f77b4', lw=2, label='Free')
ax1.plot(r_values, C_exact_r, 's--', color='#d62728', lw=2,
         label=f'Interacting (lambda={lam})')
ax1.set_xlabel('Separation r = |i-j|', fontsize=13)
ax1.set_ylabel(r'$C(r) = \langle \phi_i \phi_j \rangle$', fontsize=13)
ax1.set_title('Equal-time correlator', fontsize=13)
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)

# Right: log scale -> exponential decay shows as straight line
ax2.semilogy(r_values, abs(C_free_r), 'o-', color='#1f77b4', lw=2, label='Free')
ax2.semilogy(r_values, abs(C_exact_r), 's--', color='#d62728', lw=2,
             label=f'Interacting (lambda={lam})')
ax2.set_xlabel('Separation r = |i-j|', fontsize=13)
ax2.set_ylabel(r'$C(r)$ (log scale)', fontsize=13)
ax2.set_title('Correlator: exponential decay', fontsize=13)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('07_phi4_propagator.png', dpi=150)
print("\nSaved: 07_phi4_propagator.png")

# ================================================================
# Step 6: Extract screening mass from C(r) via cosh fit
#
# On a periodic lattice of N sites, the 2-point correlator for
# a free scalar of mass m has the form:
#   C(r) ~ A * cosh(m * (N/2 - r))
# This accounts for correlation wrapping around both directions.
#
# We fit C(r) for r >= 1 to extract the screening mass.
# ================================================================
from scipy.optimize import curve_fit

def fit_screening_mass(C, N):
    """Fit C(r) to A * cosh(m * (N/2 - r)) for r in [1, N//2-1]."""
    def cosh_model(r, A, m):
        return A * np.cosh(m * (N / 2.0 - r))
    r_data = np.arange(1, N // 2)
    try:
        popt, _ = curve_fit(cosh_model, r_data, C[r_data],
                            p0=[C[1], 0.5], maxfev=2000)
        return popt[1]  # m
    except (RuntimeError, ValueError):
        return np.nan

m_screen_free = fit_screening_mass(C_free_r, N_sites)
m_screen_exact = fit_screening_mass(C_exact_r, N_sites)

# Physical mass in free theory: lowest normal mode frequency
m_phys_free = omega_k[0]

print(f"\n--- Screening mass from correlator ---")
print(f"Free theory:")
print(f"  Lowest normal mode (physical mass): m = {m_phys_free:.4f}")
print(f"  cosh fit to C(r):                   m_screen = {m_screen_free:.4f}")
print(f"Interacting (lambda={lam}):")
print(f"  cosh fit to C(r):                   m_screen = {m_screen_exact:.4f}")

# ================================================================
# Step 7: Scan correlator over lambda
# ================================================================
lam_scan = np.array([0.0, 0.05, 0.1, 0.2, 0.5])
m_screen_scan = []

# Free GS energy
E0_free = H0_diag[0]

# Build unscaled phi_i^4 sum for re-use: H1 = lam/24 * phi4_sum
# phi4_sum = (24/lam) * H1 = sum_i phi_i^4
phi4_sum = H1 * (24.0 / lam)

print(f"\n--- Mass renormalization scan ---")
print(f"{'lambda':>8s}  {'E0':>10s}  {'<phi^2>':>10s}  {'m_screen':>10s}")
print(f"  {'0 (free)':>8s}  {E0_free:.6f}  {phi2_free.mean():.6f}  {m_phys_free:.6f}")

for lam_i in lam_scan[1:]:
    H1_i = phi4_sum * (lam_i / 24.0)
    H_i = np.diag(H0_diag) + H1_i

    E_i, psi_i = eigh(H_i)
    psi0_i = psi_i[:, 0]
    E0_i = E_i[0]

    # Correlator
    C_i = np.zeros((N_sites, N_sites))
    for ii in range(N_sites):
        phi_i_psi = phi_ops[ii] @ psi0_i
        for jj in range(N_sites):
            C_i[ii, jj] = np.dot(phi_i_psi.conj(), phi_ops[jj] @ psi0_i)

    m_screen_i = fit_screening_mass(C_i[0, :], N_sites)
    phi2_i = np.mean([C_i[i, i] for i in range(N_sites)])
    m_screen_scan.append(m_screen_i)

    print(f"  {lam_i:.4f}  {E0_i:.6f}  {phi2_i:.6f}  {m_screen_i:.6f}")

# ================================================================
# Step 8: Plot mass vs lambda
# ================================================================
fig2, ax = plt.subplots(figsize=(8, 5))

ax.axhline(y=m_phys_free, color='gray', linestyle=':', lw=1.5,
           label=f'Free mass m = {m_phys_free:.4f}')
ax.plot(lam_scan[1:], m_screen_scan, 'o-', color='#d62728', lw=2, markersize=8,
        label='Screening mass from correlator')
ax.set_xlabel(r'$\lambda$', fontsize=14)
ax.set_ylabel(r'$m_{\rm screen}$', fontsize=14)
ax.set_title(r'Mass renormalization in $\phi^4$ theory', fontsize=13)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('07_mass_renormalization.png', dpi=150)
print("\nSaved: 07_mass_renormalization.png")

# ================================================================
# Summary
# ================================================================
print()
print("=== Summary ===")
print("The equal-time 2-point correlator C_{ij} = <0|phi_i phi_j|0>")
print("is the LATTICE PROPAGATOR. In the free theory it decays")
print("exponentially with correlation length xi = 1/m.")
print("With lambda > 0, the interaction dresses the vacuum:")
print("  - <phi_i^2> changes -> mass renormalization")
print("  - The propagator shape changes")
print("This is the foundation for all Feynman diagram calculations:")
print("  D(x-y) = int d^3k/(2pi)^3 e^{ik(x-y)} / (2 omega_k)")
print("The numerical exact diagonalization goes beyond perturbation")
print("theory and captures the full non-perturbative mass shift.")
