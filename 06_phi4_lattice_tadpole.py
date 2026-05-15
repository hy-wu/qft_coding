"""
Phase 2.1: phi^4 on the Lattice — The Tadpole Sum Over Modes
=============================================================
Extend the single anharmonic oscillator to N coupled oscillators.
This IS phi^4 field theory on a 1D lattice.

Physics:
  H = sum_i [ 1/2 pi_i^2 + 1/2 m^2 phi_i^2 + 1/2 (phi_{i+1} - phi_i)^2 ]
      + lambda/4! sum_i phi_i^4

  Free part: diagonalized by Fourier transform → independent modes omega_k
  Interaction: lambda phi^4 couples modes
  1st-order GS shift: <0| H_1 |0> = lambda/4! sum_i <0| phi_i^4 |0>

  The single-site expectation <0| phi_i^4 |0> expands to a sum over modes
  via Wick contractions → the TADPOLE diagram.

Result (translation-invariant, N sites):
  Delta E0^(1) = lambda/4! * N * 3 * [1/N sum_k 1/(2 omega_k)]^2
               = lambda/8 * (1/N) * [sum_k 1/(2 omega_k)]^2

  In the continuum limit N → infinity:
  Delta E0^(1) = lambda/8 * [∫ dk/(2pi) 1/(2 omega_k)]^2

  This is the zero-temperature tadpole: the T=0 version of Exercise 01!
"""

import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# ================================================================
# Parameters
# ================================================================
m = 0.5                     # bare mass
g = 1.0                     # nearest-neighbor coupling
lam = 0.05                  # phi^4 coupling
N_sites = 4                 # number of lattice sites
N_fock = 6                  # Fock truncation per mode (0..N_fock-1)

# ================================================================
# Step 1: Diagonalize the free Hamiltonian
# H_0 = 1/2 sum_i [pi_i^2 + m^2 phi_i^2 + (phi_{i+1}-phi_i)^2]
#
# For periodic BC, the normal modes are Fourier modes.
# We'll do it numerically with exact diagonalization.
# ================================================================

# Build K matrix (phi^2 term: mass + gradient)
K = np.zeros((N_sites, N_sites))
for i in range(N_sites):
    # Nearest-neighbor coupling (periodic BC)
    K[i, i] = m**2 + 2 * g
    K[i, (i+1) % N_sites] = -g
    K[i, (i-1) % N_sites] = -g

# Diagonalize K → mode frequencies omega_k^2
omega2_k, U = eigh(K)       # U[:, k] is the k-th normal mode
omega_k = np.sqrt(omega2_k)

# Sort modes by frequency
idx_sort = np.argsort(omega_k)
omega_k = omega_k[idx_sort]
U = U[:, idx_sort]

print("=== phi^4 on 1D lattice ===")
print(f"N_sites = {N_sites},  m = {m},  g = {g},  lambda = {lam}")
print(f"\nNormal mode frequencies:")
for k in range(N_sites):
    print(f"  mode {k}: omega_{k} = {omega_k[k]:.4f}")

# ================================================================
# Step 2: Build the Hamiltonian in the product Fock basis
#
# Each mode has its own Fock space (truncated to N_fock).
# Total basis: |n_0, n_1, ..., n_{N-1}>
# Dimension: N_fock ** N_sites
# ================================================================
dim = N_fock ** N_sites
print(f"\nHilbert space dimension: {dim}")

# Map global index to Fock state tuple (multi-index)
def idx_to_state(idx):
    state = []
    for _ in range(N_sites):
        state.append(idx % N_fock)
        idx //= N_fock
    return tuple(state)

# Map Fock state tuple to global index
def state_to_idx(state):
    idx = 0
    for i, n in enumerate(state):
        idx += n * (N_fock ** i)
    return idx

# ================================================================
# Step 3: Build H_0 (diagonal in Fock basis)
# ================================================================
H0_diag = np.zeros(dim)
for idx in range(dim):
    state = idx_to_state(idx)
    E = 0.0
    for k, n in enumerate(state):
        E += omega_k[k] * (n + 0.5)
    H0_diag[idx] = E

# ================================================================
# Step 4: Build H_1 = lambda/4! sum_i phi_i^4
#
# phi_i = sum_k U[i,k] * X_k
# X_k = 1/sqrt(2 omega_k) * (a_k + a_k^dagger)
#
# H_1 = lambda/4! sum_i (sum_k U[i,k] X_k)^4
#     = lambda/4! sum_i sum_{k1,k2,k3,k4} U[i,k1]U[i,k2]U[i,k3]U[i,k4]
#       * X_{k1} X_{k2} X_{k3} X_{k4}
#
# We compute matrix elements directly in the product basis.
# ================================================================

# Precompute single-mode x matrix (truncated Fock basis)
def build_x_matrix(Nf, omega):
    """x = 1/sqrt(2*omega) (a + a^dagger) in Fock basis, truncated to Nf x Nf"""
    x = np.zeros((Nf, Nf))
    for n in range(Nf):
        if n > 0:
            x[n, n-1] = np.sqrt(n)
        if n < Nf - 1:
            x[n, n+1] = np.sqrt(n + 1)
    return x / np.sqrt(2 * omega)

# Build phi_i operator in the full Hilbert space
# by taking tensor products of single-mode operators

# Actually, let me use a more efficient approach.
# For each mode k, build the x_k operator acting on the full space.
X_ops = []
for k in range(N_sites):
    x_k_single = build_x_matrix(N_fock, omega_k[k])
    # Promote to full space by tensor product
    op_list = [np.eye(N_fock)] * N_sites
    op_list[k] = x_k_single
    X_k = op_list[0]
    for op in op_list[1:]:
        X_k = np.kron(X_k, op)
    X_ops.append(X_k)

# Build H_1 = lambda/4! sum_i (sum_k U[i,k] X_k)^4
H1 = np.zeros((dim, dim))

for i in range(N_sites):
    # phi_i operator = sum_k U[i,k] * X_k
    phi_i = sum(U[i, k] * X_ops[k] for k in range(N_sites))
    # phi_i^4
    phi_i_4 = phi_i @ phi_i @ phi_i @ phi_i
    H1 += phi_i_4

H1 *= lam / 24.0  # 1/4! = 1/24

# ================================================================
# Step 5: Exact diagonalization
# ================================================================
H_full = np.diag(H0_diag) + H1

print("\nDiagonalizing... (dim =", dim, ")")
E_exact, psi_exact = eigh(H_full)
E0_exact = E_exact[0]
psi0 = psi_exact[:, 0]

print(f"Exact GS energy:  E0 = {E0_exact:.6f}")

# ================================================================
# Step 6: 1st-order perturbation theory
# ================================================================
# GS of H_0 is |0,0,...,0>
gs_idx = state_to_idx(tuple([0] * N_sites))
E0_free = H0_diag[gs_idx]
E0_1st = E0_free + H1[gs_idx, gs_idx]

print(f"Free GS energy:   E0^(0) = {E0_free:.6f}")
print(f"1st-order PT:     E0^(1) = {E0_1st:.6f}")
print(f"  <0|H_1|0>       = {H1[gs_idx, gs_idx]:.6f}")
print(f"Error vs exact:   {abs(E0_1st - E0_exact):.2e}")

# ================================================================
# Step 7: The tadpole sum over modes
#
# Translation-invariant analytic result:
# <phi_i^4> = 3 * (1/N sum_k 1/(2 omega_k))^2
# Delta E = lam/24 * N * <phi_i^4>
# ================================================================
sum_1_over_2omega = sum(1.0 / (2 * w) for w in omega_k) / N_sites
phi2_vac = sum_1_over_2omega  # <0|phi_i^2|0> = 1/N sum_k 1/(2*omega_k)
phi4_vac = 3 * phi2_vac**2   # <0|phi_i^4|0> = 3 * (<0|phi_i^2|0>)^2

E0_tadpole_analytic = lam / 24.0 * N_sites * phi4_vac

print(f"\n--- Tadpole sum over modes ---")
print(f"<0|phi_i^2|0> = (1/N) sum_k 1/(2 omega_k) = {phi2_vac:.6f}")
print(f"<0|phi_i^4|0> = 3 * (<phi_i^2>)^2 = {phi4_vac:.6f}")
print(f"N sites * <phi_i^4> = {N_sites * phi4_vac:.6f}")
print(f"Analytic tadpole: lam/24 * N * <phi_i^4> = {E0_tadpole_analytic:.6f}")
print(f"Numeric <0|H_1|0>:  {H1[gs_idx, gs_idx]:.6f}")
print(f"Match: {abs(E0_tadpole_analytic - H1[gs_idx, gs_idx]):.2e}")

# ================================================================
# Step 8: Compare with single-mode result
# ================================================================
print(f"\n--- Single-mode comparison ---")
print(f"Single QHO:  <0|x^4|0> = 3/(4 omega^2)")
for k in range(N_sites):
    single_tad = 3.0 / (4.0 * omega_k[k]**2)
    print(f"  mode {k} (omega={omega_k[k]:.3f}):  <0|x^4|0> = {single_tad:.6f}")

# ================================================================
# Step 9: Ground state wavefunction in Fock basis
# ================================================================
print(f"\n--- Ground state composition ---")
# List the largest components
n_components = min(10, dim)
components = []
for idx in range(dim):
    amp = abs(psi0[idx])
    if amp > 0.01:
        state = idx_to_state(idx)
        components.append((amp, state))

components.sort(key=lambda x: -x[0])
for amp, state in components[:n_components]:
    print(f"  |{','.join(str(n) for n in state)}>  amplitude = {amp:.4f}")

# ================================================================
# Plot: ground state energy vs lambda
# ================================================================
lam_values = np.logspace(-2, 0, 8)
E0_exact_scan = []
E0_1st_scan = []
E0_tadpole_scan = []

for lam_i in lam_values:
    H1_i = H1 * (lam_i / lam)  # scale
    H_i = np.diag(H0_diag) + H1_i
    E_i = eigh(H_i, subset_by_index=[0, 0])[0][0]
    E0_exact_scan.append(E_i)
    E0_1st_scan.append(H0_diag[gs_idx] + H1_i[gs_idx, gs_idx])
    E0_tadpole_scan.append(H0_diag[gs_idx] + lam_i / 24.0 * N_sites * phi4_vac)

plt.figure(figsize=(8, 5))
plt.plot(lam_values, E0_exact_scan, 'o-', color='#1f77b4', lw=2, label='Exact ED')
plt.plot(lam_values, E0_1st_scan, 's--', color='#d62728', lw=2, label='1st-order PT')
plt.plot(lam_values, E0_tadpole_scan, ':^', color='#2ca02c', lw=1.5,
         label='Tadpole analytic')
plt.xscale('log')
plt.xlabel(r'$\lambda$', fontsize=14)
plt.ylabel(r'$E_0$', fontsize=14)
plt.title(rf'$\phi^4$ on {N_sites}-site lattice: GS energy ($m={m}$)', fontsize=13)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('06_phi4_lattice_E0.png', dpi=150)

# ================================================================
# Summary
# ================================================================
print()
print("=== Summary ===")
print(f"The tadpole diagram in phi^4 field theory is just the sum")
print(f"over momentum modes of the single-oscillator <0|x^4|0> = 3/(4 omega^2).")
print(f"In the continuum limit (N->inf):")
print(f"  E0_tadpole^(1) = lambda/8 * [∫ dk/(2pi) 1/(2 omega_k)]^2")
print(f"This is exactly the T=0 version of Exercise 01 (finite-T tadpole).")
print(f"The only difference is the statistical factor:")
print(f"  T=0:  1/(2 omega_k)")
print(f"  T>0:  1/(2 omega_k) + n_B(omega_k)/omega_k")
