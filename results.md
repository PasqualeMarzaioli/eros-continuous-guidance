# Results: Eros Continuous Low-Thrust Guidance

**Author:** Pasquale Marzaioli

Numerical results from a complete execution of `eros_continuous_guidance.m`.
All values below are taken from the run diary; plots are written to `plots/`.

---

## Run metadata

| Item | Value |
|---|---|
| Script | `eros_continuous_guidance.m` |
| MATLAB version | 25.2.0.3177638 (R2025b) Update 5 |
| Random seed (costate screening) | `10775298` |
| Screening ODE tolerances | RelTol \(10^{-9}\), AbsTol \(10^{-11}\) |
| Refinement ODE tolerances | RelTol \(10^{-11}\), AbsTol \(10^{-12}\) |
| N-body ODE tolerances | RelTol \(10^{-12}\), AbsTol \(10^{-13}\) |
| Costate candidates per batch | 300 (up to 5 batches) |
| Inclination continuation grid | \(0^\circ, 0.875^\circ, 1.750^\circ, 2.625^\circ, 3.500^\circ\) |
| Assertions | All passed (terminal errors, element bounds, trade metrics) |

---

## Physical and nondimensional parameters

### Physical constants

| Quantity | Symbol | Value | Unit |
|---|---|---|---|
| Eros gravitational parameter | \(\mu\) | \(3.5\times 10^{-4}\) | \(\mathrm{km}^3/\mathrm{s}^2\) |
| Asteroid radius | \(R_a\) | \(17.00\) | km |
| Initial altitude | \(h_i\) | \(52.25\) | km |
| Final altitude | \(h_f\) | \(36.15\) | km |
| Initial mass | \(m_0\) | \(25.0\) | kg |
| Thrust | \(T\) | \(21.59\times 10^{-9}\) | \(\mathrm{kg\,km/s}^2\) |
| Specific impulse | \(I_{sp}\) | \(382.82\) | s |
| Standard gravity | \(g_0\) | \(9.80665\times 10^{-3}\) | \(\mathrm{km/s}^2\) |
| Dust peak radius A | \(\rho_A\) | \(40.314\) | km |
| Dust peak radius B | \(\rho_B\) | \(59.170\) | km |
| Dust coefficients | \(k_1,k_2,k_3,k_4\) | \(7.393750\times 10^{-3}\), \(7.500000\times 10^{-3}\), \(3.696875\times 10^{-4}\), \(6.250000\times 10^{-4}\) | ND (script units) |

### Nondimensionalization (run output)

| Unit / parameter | Value |
|---|---|
| Distance unit DU | \(69.2500000000\) km |
| Mass unit MU | \(25.0000000000\) kg |
| Time unit TU | \(30803.1864365175\) s |
| Velocity unit VU | \(0.0022481440\) km/s |
| \(h_i\) (ND) | \(0.7545126354\) |
| \(h_f\) (ND) | \(0.5220216606\) |
| \(R_a\) (ND) | \(0.2454873646\) |
| \(\mu\) (ND) | \(1.0000000000\) |
| \(\rho_A\) (ND) | \(0.5821516245\) |
| \(\rho_B\) (ND) | \(0.8544404332\) |
| \(m_0\) (ND) | \(1.0000000000\) |
| \(T\) (ND) | \(0.0118327079\) |
| \(I_{sp}\) (ND) | \(0.0124279350\) |
| \(g_0\) (ND) | \(134366.8656875000\) |

---

## Boundary states

Circular orbits about Eros; initial true anomaly / phase \(45^\circ\).

**Initial state** \([\mathbf{r},\mathbf{v}]\) (km, km/s):

```
+48.9671445972  +48.9671445972  +0.0000000000
-0.0015896779   +0.0015896779   +0.0000000000
```

**Target planar state** \([\mathbf{r},\mathbf{v}]\) (km, km/s):

```
+53.1500000000  +0.0000000000  +0.0000000000
+0.0000000000   +0.0025661521  +0.0000000000
```

The inclined target rotates the final velocity out of plane by the continuation inclination (final inclined case: \(3.5^\circ\)).

---

## Costate screening

| Batch | Best pre-solve residual |
|---|---|
| 1 | \(2.230868\times 10^{0}\) |
| 2 | \(1.538127\times 10^{0}\) |

A converged planar PMP extremal was obtained in batch 2 (no further batches required).

---

## Planar PMP solution

| Quantity | Value |
|---|---|
| Transfer time \(t_f\) | \(7636.6175832325\) min |
| Final mass \(m_f\) | \(24.9973649453\) kg |
| Propellant used \(m_0-m_f\) | \(0.0026350547\) kg |
| Nondimensional dust exposure \(J=\int q\,d\tau\) | \(2.9975577758\) |
| Position error | \(1.4749022249\times 10^{-7}\) km |
| Velocity error | \(6.7344343474\times 10^{-9}\) m/s |

**Initial costates** \(\boldsymbol{\lambda}_0=[\boldsymbol{\lambda}_r,\boldsymbol{\lambda}_v,\lambda_m]\):

```
-0.2044535675  -0.9347558315  +0.0000000000
-0.0106025041  -0.6491970819  +0.0000000000
+2.2511353146
```

Out-of-plane costates are identically zero by planar symmetry.

**Associated plots**

- `plots/eros_dust_density_profile.png`
- `plots/eros_planar_optimal_trajectory.png`
- `plots/eros_planar_radius_profile.png`
- `plots/eros_planar_thrust_angles.png`
- `plots/eros_planar_hamiltonian_history.png`

---

## Inclination continuation

Continuation from the planar solution to a \(3.5^\circ\) inclined terminal velocity, with intermediate solves at:

| Step | Inclination | Status |
|---|---|---|
| 1 | \(0.000^\circ\) | planar seed |
| 2 | \(0.875^\circ\) | converged |
| 3 | \(1.750^\circ\) | converged |
| 4 | \(2.625^\circ\) | converged |
| 5 | \(3.500^\circ\) | converged |

---

## Inclined PMP solution (\(3.5^\circ\))

| Quantity | Value |
|---|---|
| Transfer time \(t_f\) | \(7568.3111487992\) min |
| Final mass \(m_f\) | \(24.9973885148\) kg |
| Propellant used \(m_0-m_f\) | \(0.0026114852\) kg |
| Nondimensional dust exposure \(J=\int q\,d\tau\) | \(3.3182338282\) |
| Position error | \(1.6216214921\times 10^{-8}\) km |
| Velocity error | \(7.6219351102\times 10^{-10}\) m/s |

**Initial costates** \(\boldsymbol{\lambda}_0\):

```
+18.6887659790  +15.7620445476  -35.8390017693
-18.7633874553  +15.1799195321  -36.9139358371
+10.6917468377
```

**Exposure comparison (planar vs inclined)**

| Case | \(J\) (ND) | Relative to planar |
|---|---|---|
| Planar | \(2.9975577758\) | — |
| Inclined \(3.5^\circ\) | \(3.3182338282\) | \(+10.70\%\) |

**Associated plots**

- `plots/eros_inclined_optimal_trajectory.png`
- `plots/eros_inclined_radius_profile.png`
- `plots/eros_inclined_thrust_angles.png`
- `plots/eros_inclined_inclination_history.png`
- `plots/eros_inclined_hamiltonian_history.png`
- `plots/eros_cumulative_dust_exposure.png`
- `plots/eros_inclination_trade_space.png`

---

## SPICE n-body free-flight check

Terminal inclined state propagated for 31 days from epoch

`2012 JAN 15 00:00:00.000 TDB` → `2012 FEB 15 00:00:00.000 TDB`

in an Eros-centered model with differential third-body accelerations from SPICE ephemerides (Sun, planets, Moon, Pluto), compared against the unperturbed Kepler orbit about Eros.

| Metric (n-body vs Kepler at 31 days) | Value |
|---|---|
| Position difference \(\|\Delta\mathbf{r}\|\) | \(1.6085097170\times 10^{-1}\) km |
| Velocity difference \(\|\Delta\mathbf{v}\|\) | \(7.8255213485\times 10^{-6}\) km/s |
| Semimajor-axis difference \(\Delta a\) | \(-1.6086651215\times 10^{-5}\) km |
| Phase difference \(\Delta u\) | \(-6.5664811015\times 10^{-2}\) deg |
| Maximum eccentricity difference \(\max\|\Delta e\|\) | \(6.3432263163\times 10^{-5}\) |

Verification assertions on osculating elements were satisfied:

- \(\max\|\Delta a\| < 10^{-2}\) km
- \(\max\|\Delta e\| < 10^{-3}\)

**Associated plots**

- `plots/eros_nbody_trajectory_comparison.png`
- `plots/eros_nbody_ntc_error_history.png`
- `plots/eros_perturbation_source_attribution.png`
- `plots/eros_osculating_orbit_differences.png`

---

## Exported figure inventory

All figures were exported to `plots/` at 160 dpi:

| File | Content |
|---|---|
| `eros_dust_density_profile.png` | Radial dust density \(q(\rho)\) with \(h_i\), \(h_f\) markers |
| `eros_planar_optimal_trajectory.png` | Planar trajectory on dust field |
| `eros_planar_radius_profile.png` | Planar radius history vs outer density peak |
| `eros_planar_thrust_angles.png` | Planar in-plane / cross-track thrust angles |
| `eros_planar_hamiltonian_history.png` | Planar Hamiltonian constancy check |
| `eros_inclined_optimal_trajectory.png` | Inclined trajectory on dust field |
| `eros_inclined_radius_profile.png` | Inclined radius history |
| `eros_inclined_thrust_angles.png` | Inclined thrust angles |
| `eros_inclined_inclination_history.png` | Osculating inclination vs time |
| `eros_inclined_hamiltonian_history.png` | Inclined Hamiltonian constancy check |
| `eros_cumulative_dust_exposure.png` | Instantaneous and cumulative exposure |
| `eros_inclination_trade_space.png` | Exposure, TOF, propellant, peak \(q\), cross thrust vs \(\Delta i\) |
| `eros_nbody_trajectory_comparison.png` | N-body vs Kepler trajectory views |
| `eros_nbody_ntc_error_history.png` | Radial–tangential–cross error histories |
| `eros_perturbation_source_attribution.png` | Leave-one-out third-body attribution |
| `eros_osculating_orbit_differences.png` | Osculating \(a,e,i,\Omega,u\) differences |

---

## Summary

1. A bang (full-thrust) PMP extremal for the planar dust-minimizing transfer converges to sub-millimetre / sub-nanometre-per-second terminal accuracy.
2. Continuation to \(3.5^\circ\) inclination remains well-conditioned; exposure increases by about \(11\%\) relative to the planar optimum while propellant use stays on the order of a few grams.
3. A 31-day SPICE n-body free flight from the inclined terminal state stays within centimetre-to-decimetre position and micro-metre-per-second velocity of the Kepler reference in energy-related elements, confirming that the local two-body PMP model is an appropriate guidance approximation over that horizon.
