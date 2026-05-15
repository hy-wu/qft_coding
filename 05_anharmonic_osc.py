"""
Phase 2.0: Anharmonic Oscillator — Perturbation Theory vs Exact Diagonalization
==============================================================================

Physics:
  H = H_0 + H_1
  H_0 = 1/2 p^2 + 1/2 omega^2 x^2     (free QHO)
  H_1 = lambda x^4                       (anharmonic perturbation)

We compute the ground state energy shift in two ways:
  1. Rayleigh-Schrodinger perturbation theory (1st and 2nd order)
  2. Exact diagonalization in a truncated Fock basis

The 1st-order shift Delta E_0^(1) = <0| lambda x^4 |0>
is the ZERO-TEMPERATURE TADPOLE diagram — the T=0 version
of our Exercise 01 finite-T tadpole!

This is the foundation: in phi^4 field theory, each site's
lambda phi(x)^4 couples modes → the tadpole is just this
anharmonic ground state shift, summed over all momenta.
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal
import matplotlib.pyplot as plt

# ================================================================
# Parameters
# ================================================================
omega = 1.0                     # harmonic frequency
lam = 0.1                       # anharmonic coupling
N_max = 30                      # Fock space truncation (|0> to |N_max-1>)

# ================================================================
# Free QHO: H_0 matrix in the Fock basis (number states)
# H_0 |n> = omega (n + 1/2) |n>
# ================================================================
E0_free = np.array([omega * (n + 0.5) for n in range(N_max)])

# ================================================================
# x operator in the Fock basis
#
# x = 1/sqrt(2 omega) (a + a^dagger)
# <n| x |m> = 1/sqrt(2 omega) * (sqrt(m) delta_{n,m-1} + sqrt(n) delta_{n-1,m})
# ================================================================
x_matrix = np.zeros((N_max, N_max))
for n in range(N_max):
    if n > 0:
        x_matrix[n, n-1] = np.sqrt(n)       # a |n> = sqrt(n) |n-1>
    if n < N_max - 1:
        x_matrix[n, n+1] = np.sqrt(n + 1)   # a^dagger |n> = sqrt(n+1) |n+1>
x_matrix /= np.sqrt(2 * omega)

# H_1 = lambda x^4
x4_matrix = x_matrix @ x_matrix @ x_matrix @ x_matrix
H1_matrix = lam * x4_matrix

# Full Hamiltonian
H_matrix = np.diag(E0_free) + H1_matrix

# ================================================================
# Exact diagonalization
# ================================================================
E_exact, states_exact = np.linalg.eigh(H_matrix)

# Ground state
E0_exact = E_exact[0]
psi0_exact = states_exact[:, 0]

print("=== Anharmonic Oscillator: lambda =", lam, "===")
print(f"Free ground state energy:  E0_free = {E0_free[0]:.6f}")
print(f"Exact ground state energy: E0_exact = {E0_exact:.6f}")

# ================================================================
# 1st-order perturbation theory
# Delta E0^(1) = <0| H_1 |0> = lambda * <0| x^4 |0>
# ================================================================
# <0| x^4 |0> = <0| x^2 |0>^2 + contractions
# Using Wick's theorem (or brute force):
# x^4 = 1/(4 omega^2) (a + a^dagger)^4
# <0| x^4 |0> = 3/(4 omega^2) = 3/(4 omega^2)
# ... let me just compute it numerically:

E0_1st_order = E0_free[0] + H1_matrix[0, 0]
print(f"\n1st-order PT:")
print(f"  <0| H_1 |0> = {H1_matrix[0, 0]:.6f}")
print(f"  E0^(1) = {E0_1st_order:.6f}")
print(f"  error vs exact: {abs(E0_1st_order - E0_exact):.6e}")

# ================================================================
# 2nd-order perturbation theory
# Delta E0^(2) = sum_{m != 0} |<m| H_1 |0>|^2 / (E0 - Em)
# ================================================================
E0_2nd_order = E0_1st_order
for m in range(1, N_max):
    denom = E0_free[0] - E0_free[m]
    if abs(denom) > 1e-15:
        E0_2nd_order += abs(H1_matrix[m, 0])**2 / denom

print(f"\n2nd-order PT:")
print(f"  E0^(2) = {E0_2nd_order:.6f}")
print(f"  error vs exact: {abs(E0_2nd_order - E0_exact):.6e}")

# ================================================================
# Scan over coupling strength
# ================================================================
lam_values = np.logspace(-2, 0.5, 12)
E0_exact_scan = []
E0_1st_scan = []
E0_2nd_scan = []

print(f"\n--- Scanning lambda ---")
for lam_i in lam_values:
    H1_i = lam_i * x4_matrix
    H_i = np.diag(E0_free) + H1_i
    E_i = np.linalg.eigvalsh(H_i)

    E0_exact_scan.append(E_i[0])

    # 1st order
    E0_1st = E0_free[0] + H1_i[0, 0]
    E0_1st_scan.append(E0_1st)

    # 2nd order
    E0_2nd = E0_1st
    for m in range(1, N_max):
        denom = E0_free[0] - E0_free[m]
        if abs(denom) > 1e-15:
            E0_2nd += abs(H1_i[m, 0])**2 / denom
    E0_2nd_scan.append(E0_2nd)

    print(f"lam={lam_i:.4f}:  exact={E_i[0]:.6f},  1st={E0_1st:.6f},  2nd={E0_2nd:.6f}")

# ================================================================
# Plot: ground state energy vs lambda
# ================================================================
plt.figure(figsize=(8, 5))

plt.plot(lam_values, E0_exact_scan, 'o-', color='#1f77b4', linewidth=2,
         label='Exact ED', markersize=6)
plt.plot(lam_values, E0_1st_scan, 's--', color='#d62728', linewidth=2,
         label='1st-order PT', markersize=6)
plt.plot(lam_values, E0_2nd_scan, '^-.', color='#2ca02c', linewidth=2,
         label='2nd-order PT', markersize=6)

plt.xscale('log')
plt.xlabel(r'$\lambda$', fontsize=14)
plt.ylabel(r'$E_0$', fontsize=14)
plt.title(r'Anharmonic oscillator: $H = \frac{1}{2}(p^2+\omega^2 x^2) + \lambda x^4$',
          fontsize=12)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('05_anh_osc_E0.png', dpi=150)

# ================================================================
# Plot: ground state wavefunction in Fock basis
# ================================================================
plt.figure(figsize=(8, 5))

n_states = np.arange(15)
psi0_overlap = np.abs(states_exact[:, 0][:15])

# Compare free vs exact ground state
psi0_free = np.zeros(15)
psi0_free[0] = 1.0  # free ground state is |0>

x_plot = n_states
width = 0.35
plt.bar(x_plot - width/2, psi0_free, width, alpha=0.6, color='gray',
        label='Free GS |0>')
plt.bar(x_plot + width/2, psi0_overlap, width, alpha=0.8, color='#1f77b4',
        label=f'Exact GS (lambda={lam})')

plt.xlabel(r'$n$ (Fock state)', fontsize=14)
plt.ylabel(r'$|\langle n | \psi_0 \rangle|$', fontsize=14)
plt.title('Ground state composition in Fock basis', fontsize=13)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig('05_anh_osc_wavefunction.png', dpi=150)

# ================================================================
# The tadpole connection
# ================================================================
print()
print("--- Tadpole connection ---")
# <0| x^4 |0> = 3/(4 omega^2)
tadpole_analytic = 3.0 / (4.0 * omega**2)
tadpole_numeric = x4_matrix[0, 0]
print(f"Analytic <0|x^4|0> = {tadpole_analytic:.6f}")
print(f"Numeric  <0|x^4|0> = {tadpole_numeric:.6f}")
print(f"1st-order shift = lambda * <0|x^4|0> = {lam * tadpole_numeric:.6f}")
print()
print("In phi^4 field theory, the ground state energy shift")
print("from the tadpole diagram is exactly this: summed over all modes.")
print("At finite T (Exercise 01), the thermal tadpole replaces")
print("the vacuum expectation with the thermal average: n_B(E_k)/E_k vs 1/(2E_k).")
