"""
Phase 1.1: Coupled Quantum Harmonic Oscillator Chain
====================================================
A discrete version of a free scalar field.

We take N QHOs on a 1D lattice with nearest-neighbor coupling
and diagonalize the Hamiltonian exactly.

The dispersion relation:
  omega_k^2 = omega_0^2 + 4g * sin^2(ka/2)

In the continuum limit (small k):
  omega_k -> sqrt(omega_0^2 + k^2)

This IS the relativistic dispersion of a free scalar field!
"""

import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# ================================================================
# Parameters
# ================================================================
N = 32                  # number of sites
omega_0 = 0.5           # on-site frequency
g = 1.0                 # nearest-neighbor coupling
a = 1.0                 # lattice spacing

# We work in dimensionless units (mass m = 1, hbar = 1)

# ================================================================
# Hamiltonian in matrix form
#
# H = sum_i [ 1/2 p_i^2 + 1/2 omega_0^2 x_i^2 ]
#     + g/2 sum_{<ij>} (x_i - x_j)^2
#
# In matrix form (explicitly writing the quadratic forms):
#   H = 1/2 p^T I p + 1/2 x^T K x
#
# where K_ij = omega_0^2 * delta_ij + g * (2 delta_ij - delta_{i,j+1} - delta_{i+1,j})
# ================================================================

# Build the "spring-constant" matrix K
K = np.zeros((N, N))
for i in range(N):
    K[i, i] = omega_0**2 + 2 * g        # on-site + two springs (left and right)
    if i > 0:
        K[i, i-1] = -g                   # coupling to left neighbor
        K[i-1, i] = -g                   # symmetric
# For a chain (not ring), the first and last sites have only one neighbor
# Our loop already handles this because if i=0, no left neighbor, etc.
# Actually, let me fix: the above counts 2g for all sites but edges should have only 1g
# Let me rebuild correctly:

K = np.zeros((N, N))
for i in range(N):
    # Coupling to left neighbor
    if i > 0:
        K[i, i-1] = -g
    # Coupling to right neighbor
    if i < N - 1:
        K[i, i+1] = -g
    # On-site: each site couples to itself with:
    # omega_0^2 + (number of neighbors) * g
    n_neighbors = 2
    if i == 0 or i == N - 1:
        n_neighbors = 1
    K[i, i] = omega_0**2 + n_neighbors * g

print("K matrix (first 5x5):")
print(K[:5, :5])

# ================================================================
# Diagonalize: K's eigenvalues are omega_k^2
#
# H = 1/2 sum_k (|p_k|^2 + omega_k^2 |x_k|^2)
#   = sum_k omega_k (a_k^dagger a_k + 1/2)
# ================================================================

eigenvalues, eigenvectors = eigh(K)
omega_k = np.sqrt(eigenvalues)

# Sort by k (momentum)
# eigenvalues already sorted by eigh (ascending)
# Wavevectors for 1D chain with open boundary
k_values = np.linspace(0, np.pi, N)   # approximate

# Analytical dispersion for infinite chain:
# omega^2(k) = omega_0^2 + 4g * sin^2(k/2)
k_analytical = np.linspace(0, np.pi, 200)
omega_analytical = np.sqrt(omega_0**2 + 4 * g * np.sin(k_analytical / 2)**2)

print(f"\nLowest eigenfrequency: omega_0 = {omega_k[0]:.4f}")
print(f"Highest eigenfrequency: omega_max = {omega_k[-1]:.4f}")
print(f"Continuum limit: sqrt(omega_0^2 + k^2)")

# ================================================================
# Create the normal mode operators (annihilation operators)
#
# For each mode k (eigenvector v_k):
#   a_k = sqrt(omega_k/2) (x_k + i/omega_k p_k)
# where
#   x_k = v_k^T x  (projection of x onto mode k)
#   p_k = v_k^T p
#
# In the original basis:
#   a_k = sqrt(omega_k/2) v_k^T x + i/(sqrt(2 omega_k)) v_k^T p
#
# The matrix representation of a_k in the position basis is...
# Well, we can't easily show the operator matrix here (it's NxN infinite).
# But we can show the normal mode coordinate - how each site contributes.
# ================================================================

# Plot the first few normal mode shapes (eigenvectors)
mode_indices = [0, 1, 2, 3, 7, 15, 31]
plt.figure(figsize=(10, 6))
for idx in mode_indices[:4]:   # first 4 modes
    if idx < N:
        mode_shape = eigenvectors[:, idx]
        plt.plot(range(N), mode_shape, 'o-', label=f'mode {idx}, ω = {omega_k[idx]:.3f}')

plt.xlabel('Site index', fontsize=13)
plt.ylabel('Mode amplitude', fontsize=13)
plt.title(f'Normal mode shapes (N={N} coupled QHOs)', fontsize=13)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('04_normal_modes.png', dpi=150)

# ================================================================
# Plot the dispersion relation
# ================================================================
plt.figure(figsize=(10, 6))

# Numerical eigenfrequencies
plt.plot(k_values, omega_k, 'o', color='#1f77b4', markersize=5,
         label=f'Numerical (N={N})', zorder=3)

# Analytical dispersion for infinite chain
plt.plot(k_analytical, omega_analytical, '-', color='#d62728', linewidth=2,
         label=r'$\omega(k) = \sqrt{\omega_0^2 + 4g\sin^2(k/2)}$')

# Continuum (relativistic) limit: omega = sqrt(omega_0^2 + k^2)
k_cont = np.linspace(0, 1.5, 100)   # only valid for small k
omega_cont = np.sqrt(omega_0**2 + k_cont**2)
plt.plot(k_cont, omega_cont, '--', color='#2ca02c', linewidth=2,
         label=r'Continuum: $\sqrt{\omega_0^2 + k^2}$')

plt.axhline(y=omega_0, color='gray', linestyle=':', alpha=0.5)
plt.axvline(x=np.pi, color='gray', linestyle=':', alpha=0.5)

plt.xlabel(r'Wavenumber $k$', fontsize=13)
plt.ylabel(r'Frequency $\omega(k)$', fontsize=13)
plt.title(f'Dispersion relation: coupled QHO chain', fontsize=13)
plt.legend(fontsize=11)
plt.xlim(0, np.pi + 0.3)
plt.ylim(0, omega_k[-1] * 1.1)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('04_dispersion.png', dpi=150)

# ================================================================
# Check: continuum limit
# small k -> omega_k^2 = omega_0^2 + k^2 + O(k^4)
# ================================================================
print("\n--- Small-k check (acoustic/relativistic regime) ---")
small_k_idx = 1  # skip the k=0 mode
if small_k_idx < N:
    k_val = k_values[small_k_idx]
    omega_num = omega_k[small_k_idx]
    omega_rel = np.sqrt(omega_0**2 + k_val**2)
    print(f"k = {k_val:.6f}")
    print(f"omega_num    = {omega_num:.6f}")
    print(f"omega_rel    = {omega_rel:.6f}")
    print(f"relative diff = {abs(omega_num - omega_rel)/omega_num:.2e}")

print(f"\nNumber of sites N = {N}")
print(f"Cuto energy = {omega_k[-1]:.4f} (regulated by lattice)")

# ================================================================
# Summary
# ================================================================
print()
print("--- Summary ---")
print(f"Coupled QHO chain: non-relativistic massive particles on a lattice")
print(f"In the continuum limit (ka << 1): omega_k ~ sqrt(m^2 + k^2)")
print(f"This is the relativistic dispersion of a free scalar field!")
print(f"The lattice provides a natural UV cutof at k ~ pi/a")
