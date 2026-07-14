<!--
results.md — Records the verified numerical run and audit conclusions.
Values come from the complete deterministic MATLAB execution described below.
Author: Pasquale Marzaioli
-->

# Verified Results and Audit Record

**Author:** Pasquale Marzaioli

## Run metadata

| Item | Value |
|---|---|
| Driver | `eros_continuous_guidance.m` |
| MATLAB | R2025b Update 5 |
| Random seed | `10775298` |
| Random batches | 2, both fully evaluated |
| Candidates per batch | 300 |
| Refined candidates per batch | 12 |
| Guidance ODE tolerance | RelTol $10^{-11}$, AbsTol $10^{-12}$ |
| Screening ODE tolerance | RelTol $10^{-9}$, AbsTol $10^{-11}$ |
| SPICE comparison tolerance | RelTol $10^{-12}$, AbsTol $10^{-13}$ |
| Physical-scaling continuation | 33 nodes, 32 successful increments |
| Inclination continuation | $0,0.875,1.750,2.625,3.500^\circ$ |
| Exported figures | 16 |
| Full-run exit status | 0 |

## Physical parameters

| Quantity | Value | Unit / meaning |
|---|---:|---|
| Eros GM | $4.4630000000\times10^{-4}$ | $\mathrm{km^3/s^2}$ |
| Reference radius | 17.0000 | km, spherical coordinate convention |
| Initial reference radius | 69.2500 | km |
| Final reference radius | 53.1500 | km |
| Initial mass | 25.0000 | kg |
| Fixed thrust | 27.530334286 | $\mu\mathrm{N}$ |
| Specific impulse | 382.82 | s |
| Standard gravity | $9.80665\times10^{-3}$ | $\mathrm{km/s^2}$ |
| Mass flow | $-7.333244924\times10^{-9}$ | kg/s |
| Synthetic peak radius A | 40.314 | km |
| Synthetic peak radius B | 59.170 | km |

The thrust is scaled from the legacy pair
$(3.5\times10^{-4}\ \mathrm{km^3/s^2},21.59\ \mu\mathrm{N})$ so the original
dimensionless thrust-to-gravity ratio is preserved after correcting Eros GM.

The radial cost coefficients are

$$
(k_1,k_2,k_3,k_4)=
(7.393750\times10^{-3},7.500000\times10^{-3},
3.696875\times10^{-4},6.250000\times10^{-4}).
$$

They define a synthetic dimensionless field, not measured dust density.

## Nondimensional parameters

| Parameter | Value |
|---|---:|
| Distance unit $DU$ | 69.2500000000 km |
| Mass unit $MU$ | 25.0000000000 kg |
| Time unit $TU$ | 27278.2322807662 s |
| Velocity unit $VU$ | 0.0025386542 km/s |
| $h_i/DU$ | 0.7545126354 |
| $h_f/DU$ | 0.5220216606 |
| $R_\mathrm{ref}/DU$ | 0.2454873646 |
| $\rho_A/DU$ | 0.5821516245 |
| $\rho_B/DU$ | 0.8544404332 |
| Thrust | 0.0118327079 |
| Specific impulse | 0.0140339006 |
| Standard gravity | 105373.9704024759 |

## Fixed boundary states

The initial state uses a $45^\circ$ phase. Values are
$[\mathbf r,\mathbf v]$ in km and km/s:

~~~text
Initial
+48.9671445972  +48.9671445972  +0.0000000000
-0.0017950996   +0.0017950996   +0.0000000000

Planar target
+53.1500000000  +0.0000000000  +0.0000000000
+0.0000000000   +0.0028977560  +0.0000000000
~~~

For inclination continuation, the target velocity is rotated about the target
radial direction while its circular-speed magnitude is retained. These are fixed
Cartesian endpoint states, not orbit-manifold constraints.

## Deterministic screening and continuation

| Random batch | Best pre-solve residual |
|---:|---:|
| 1 | 2.230868 |
| 2 | 1.538127 |

Both batches and their twelve selected refinements were considered before
retaining the minimum-cost converged anchor. Batch 2 supplied that anchor with
cost 2.9975577760 under the legacy scaling. Continuation then converged at every
physical-GM node from
$3.5\times10^{-4}$ through $4.463\times10^{-4}\ \mathrm{km^3/s^2}$.

## Planar PMP extremal

| Quantity | Value |
|---|---:|
| Transfer time | 6762.7320990774 min |
| Transfer time | 4.696341735 days |
| Final mass | 24.9970244337 kg |
| Propellant | 0.0029755663 kg |
| Dimensionless cost $J$ | 2.9975385712 |
| Terminal position error | $7.2639205835\times10^{-8}$ km |
| Terminal velocity error | $3.7083348939\times10^{-9}$ m/s |
| Maximum $|H|$ | $1.6465857521\times10^{-10}$ |

Initial costates:

~~~text
-0.2044690692  -0.9347893251  +0.0000000000
-0.0105868865  -0.6492342236  +0.0000000000
+2.2511435035
~~~

The two zero out-of-plane costates follow planar symmetry.

## Inclined PMP extremal

| Quantity | Value |
|---|---:|
| Terminal inclination | $3.500^\circ$ |
| Transfer time | 6702.2372384658 min |
| Transfer time | 4.654331416 days |
| Final mass | 24.9970510512 kg |
| Propellant | 0.0029489488 kg |
| Dimensionless cost $J$ | 3.3181557629 |
| Cost change from planar branch | $+10.69601555\%$ |
| Terminal position error | $1.1059135614\times10^{-8}$ km |
| Terminal velocity error | $5.9653280650\times10^{-10}$ m/s |
| Maximum $|H|$ | $4.6068064459\times10^{-10}$ |

Initial costates:

~~~text
+18.6807354380  +15.7533539005  -35.8261663157
-18.7550746158  +15.1712280091  -36.9006564934
+10.6884316791
~~~

Every intermediate inclination solve converged.

## SPICE third-body sensitivity

The inclined terminal state is propagated in `ECLIPJ2000` from
`2012 JAN 15 00:00:00 TDB` through `2012 FEB 15 00:00:00 TDB`. The comparison
adds differential point-mass gravity from Mercury, Venus, Earth, Moon, Mars,
Jupiter, Saturn, Uranus, Neptune, Pluto, and the Sun to point-mass Eros gravity.

| Third-body minus Kepler metric | Verified value |
|---|---:|
| Final position norm | $1.4284097379\times10^{-1}$ km |
| Final velocity norm | $7.7860890455\times10^{-6}$ km/s |
| Final semimajor-axis difference | $-3.0575615341\times10^{-5}$ km |
| Final argument-of-latitude difference | $-5.9235948645\times10^{-2}$ deg |
| Maximum eccentricity difference | $4.9712651976\times10^{-5}$ |

Equivalent headline units are 142.841 m position and 7.786 mm/s velocity. The
semimajor-axis endpoint difference is -0.0306 m; its osculating history contains
larger periodic excursions shown in the plot.

The Sun-only history nearly reproduces the secular tangential difference, while
the all-except-Sun history remains much smaller. This conclusion applies only to
the included third-body terms at the chosen epoch and inertial orientation.

## Verification evidence

| Check | Result |
|---|---|
| MATLAB Code Analyzer over driver, functions, and test | `CHECKCODE_TOTAL=0` |
| Canonical 14-by-14 Jacobian vs central differences | Passed |
| Hamiltonian gradient vs central differences | Passed |
| Full shooting Jacobian vs five-point finite differences | Relative Frobenius error $7.7350\times10^{-7}$ |
| Shooting residual reconstructed from rounded published values | $8.3774\times10^{-9}$ |
| Identical-history NTC transform | Exact zero |
| Eight-condition shooting residual | $<10^{-8}$ for accepted solves |
| Planar and inclined endpoint assertions | Passed |
| Planar and inclined Hamiltonian assertions | Passed |
| Guidance GM vs loaded `gm_de440.tpc` | Passed |
| Kepler semimajor-axis conservation check | Passed |
| SPICE element finiteness and bounds | Passed |
| Full deterministic run | Exit status 0 |
| Visual inspection of all 16 PNGs | Passed |

The derivative test is intentionally small and dependency-free. It checks the
highest-risk analytical derivatives and the corrected NTC function signature; it
does not replace the end-to-end study.

## Figure inventory

| File | Verified interpretation |
|---|---|
| `eros_dust_density_profile.png` | Dimensionless synthetic radial cost |
| `eros_planar_optimal_trajectory.png` | Planar stationary extremal on the cost field |
| `eros_planar_radius_profile.png` | Planar radius and outer cost peak |
| `eros_planar_thrust_angles.png` | Continuous planar steering angles |
| `eros_planar_hamiltonian_history.png` | Planar Hamiltonian numerical history |
| `eros_inclined_optimal_trajectory.png` | Inclined stationary extremal projection |
| `eros_inclined_radius_profile.png` | Inclined radius and outer cost peak |
| `eros_inclined_thrust_angles.png` | Inclined steering angles |
| `eros_inclined_inclination_history.png` | Osculating inclination history |
| `eros_inclined_hamiltonian_history.png` | Inclined Hamiltonian numerical history |
| `eros_cumulative_dust_exposure.png` | Instantaneous and cumulative synthetic cost |
| `eros_inclination_trade_space.png` | Metrics along the continued local branch |
| `eros_nbody_trajectory_comparison.png` | Third-body and Kepler paths at plotting scale |
| `eros_nbody_ntc_error_history.png` | Absolute rotating-frame and inertial differences |
| `eros_perturbation_source_attribution.png` | Signed NTC and Sun/non-Sun comparison |
| `eros_osculating_orbit_differences.png` | Differences in $a,e,i,\Omega,u$ |

The historical `optimal_trajectory` filenames are retained for stable links; the
figures and documentation correctly describe them as stationary extremals.

## Audit disposition

| Finding | Correction |
|---|---|
| Eros GM inconsistent with bundled kernel and published measurement | Set to $4.463\times10^{-4}\ \mathrm{km^3/s^2}$; runtime kernel assertion added |
| Corrected GM reduced control authority when thrust was left unchanged | Thrust scaled to preserve the original dimensionless ratio; assumption documented |
| Continuation used the final thrust instead of each node's thrust | Corrected both nondimensionalization references |
| Invalid `ntcDifferences` declaration | Function signature repaired and regression added |
| First feasible random batch could terminate the search | Every configured batch is evaluated; lowest-cost converged candidate retained |
| Synthetic field described with physical density/exposure units | Relabeled as dimensionless residence-time cost proxy |
| Fixed-throttle control described as bang-bang | Corrected to direction-only full-throttle PMP control |
| PMP necessary conditions presented as proof of optimality | Results classified as local stationary candidates |
| Fixed terminal state described as a terminal orbit manifold | Endpoint scope corrected |
| Third-body comparison presented as point-mass validation | Reclassified as sensitivity analysis with explicit omissions |
| Position and velocity unit conversions overstated accuracy | Corrected to metres and millimetres per second |
| Plot labels, phase wrapping, legends, and blank attribution tile | Corrected and visually inspected |
| MICE described as vendored | Documented as an ignored local official dependency |
| No derivative regression | Added `tests/verify_derivatives.m` |

## Scientific interpretation

The software correctly solves the stated simplified PMP problem to tight numerical
tolerances. The result is useful as a reproducible algorithm demonstration and as
an initializer. It is not yet a physically complete near-Eros mission design:
Eros gravity harmonics/polyhedral gravity, rotation, solar radiation pressure,
shape clearance, uncertainty, and closed-loop control remain outside the force and
guidance models.
