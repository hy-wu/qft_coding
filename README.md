# qft_coding — Finite Temperature QFT with Numerical Computing

Learning QFT through computation. From Matsubara sums to spectral functions.

Current exercises:

## 01 — Tadpole & Thermal Mass

**Content:** One-loop tadpole self-energy in finite-temperature $\phi^4$ theory.

- Matsubara sum performed analytically → numerical momentum integration
- Verification of high-temperature limit $m_{\rm th}^2 = \lambda T^2 / 24$
- Plot: $m_{\rm th}^2(T)$ across $T/m_0 \in [0.3, 30]$

**Run:** `python 01_tadpole_thermal_mass.py`

## 02 — Spectral Function & Thermal Broadening

**Content:** Scalar spectral function $\rho(\omega,p)$ at finite $T$ with Breit-Wigner form.

- Thermal width $\Gamma \sim \lambda^2 T$ from $2\to2$ scattering
- Three plots: $T$-scan, $p$-scan, peak position/width vs $T$
- Physics analogue: quarkonium ($J/\psi$, $\Upsilon$) melting in QGP

**Run:** `python 02_spectral_function.py`

## Roadmap

| # | Topic | Status |
|---|-------|--------|
| 1 | Tadpole: thermal mass | ✅ done |
| 2 | Spectral function | ✅ done |
| 3 | HTL resummation / heavy quark diffusion | ⏳ next |

## Dependencies

- numpy, scipy, matplotlib

## References

- Kapusta & Gale, *Finite-Temperature Field Theory*
- Laine & Vuorinen, *Basics of Thermal Field Theory*
- Le Bellac, *Thermal Field Theory*
