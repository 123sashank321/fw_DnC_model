# fw_DnC_model

Longitudinal dynamics and control for a fixed-wing UAV, in pure MATLAB.

A complete implementation of the design chain in Beard & McLain, *Small Unmanned
Aircraft: Theory and Practice*, chapters 3 to 6, reduced to the **longitudinal
plane**: nonlinear plant, trim, linearisation, transfer functions, reduced-order
modes, and the successive-loop-closure autopilot.

Two airframes, switchable from one line: the book's **Aerosonde** and a **canard**
configuration.

Coursework for IITB AE695 (State Space Analysis and Design).

## Quick start

```matlab
run_lon        % the whole chain: printout, modal table, open- and closed-loop plots
check_lon      % 37 assertions against the book and against itself; throws on failure
```

No toolboxes required. Base MATLAB only — trim uses `fzero` rather than `fmincon`,
and transfer functions are returned as `{num, den}` coefficient vectors rather
than `tf` objects. If you want LTI objects, `tf(coef.T_theta_delta_e{:})` is
yours to call.

Switch airframe by editing one line at the top of `run_lon.m`:

```matlab
airframe = 'aerosonde';          % or 'canard'
```

## State and input

```
x = [u; w; q; theta; h]        m/s, m/s, rad/s, rad, m
u = [delta_e; delta_t]         rad, 0..1
```

Five states, matching book equation 5.50. Altitude `h` is integrated directly
rather than `p_d`, which means `A_lon` comes out in the book's `[u w q theta h]`
ordering with no row or column sign flip anywhere.

## Files

| File | Chapter | Purpose |
|---|---|---|
| `lon_params.m` | — | Both airframes in one canonical struct. The only place numbers live. |
| `f_lon.m` | 3 + 4 | Nonlinear EOM, forces and moments, both thrust models. |
| `compute_lon_trim.m` | 5.3 / F.2 | Trim at a given airspeed and flight path angle. |
| `compute_lon_ss.m` | 5.5 / F.4 | Finite-difference Jacobian → `A_lon`, `B_lon`. |
| `compute_lon_tf.m` | 5.4 | The five transfer functions and their `a`-coefficients. |
| `lon_modes.m` | 5.6 | Short-period and phugoid, exact vs. the approximations. |
| `lon_gains.m` | 6.4 | Gain design by successive loop closure. |
| `lon_autopilot.m` | 6.4 | Four loops plus the altitude state machine. |
| `lon_autopilot_reset.m` | 6.4 | Fresh integrator state. |
| `run_lon.m` | — | Driver. |
| `check_lon.m` | — | The self-check. |

`run_lon` writes `lon_coef.mat` (gitignored, regenerable). It is a cache for
working across sessions, not a hidden dependency: `lon_gains` takes `coef` as an
argument and asserts the airframe tag matches, so a stale file fails loudly
instead of quietly designing gains for the wrong aircraft.

## The transfer functions

| | Book eq. | Form |
|---|---|---|
| pitch from elevator | 5.29 | `a_theta3 / (s^2 + a_theta1 s + a_theta2)` |
| altitude from pitch | 5.31 | `Va / s` |
| altitude from airspeed | 5.32 | `theta / s` |
| airspeed from throttle | 5.36 | `a_V2 / (s + a_V1)` |
| airspeed from pitch | 5.36 | `-a_V3 / (s + a_V1)` |

## Why trim is a root-find, not an optimisation

Book Appendix F.2 shows the wings-level longitudinal problem is fully
**triangular**, so nothing needs solving simultaneously:

1. `alpha` fixes the kinematics: `theta = alpha + gamma`, `q = 0`,
   `u = Va cos(alpha)`, `w = Va sin(alpha)`.
2. `delta_e` is then **closed form** from `q_dot = 0`:
   `delta_e = -(C_m_0 + C_m_alpha*alpha) / C_m_delta_e`.
3. `alpha` is a scalar root-find on `w_dot = 0`. Thrust acts along body x and so
   never enters `f_z`, which is what makes this step independent of throttle.
4. `delta_t` is a scalar root-find on `u_dot = 0`, done last.

`h_dot` needs no equation: `u sin(theta) - w cos(theta) = Va sin(gamma)` holds
identically once `theta = alpha + gamma`, so the climb rate comes out for free.

The book's own companion code uses `fmincon` on a 17-element vector, because it
trims all twelve states at once. For the longitudinal plane that is a constrained
nonlinear program solving a triangular problem. `fzero` reaches a residual of
**~2e-15** here, against the ~1e-6 you would tune `fmincon` to, and it drops the
Optimization Toolbox dependency.

## Validation

`check_lon` runs 37 assertions in four groups. All pass.

**1. Trim is a true equilibrium.** `|f(x*,u*)| < 1e-9` for both airframes at
`gamma = 0` and `gamma = 0.1`; actual residuals are ~2e-15. The climbing case is
the one that catches a gravity or `theta` sign error, because such an error
cancels at `gamma = 0`.

**2. Two independent derivations agree.** The numerical Jacobian against the
closed-form section 5.4 coefficients: `A(3,3) = -a_theta1` and
`B(3,1) = a_theta3` hold exactly, for both airframes. This is a stronger
statement than matching a stored number, and unlike group 3 it also covers the
canard, for which no published values exist.

**3. Against published reference values** from the book's companion repository
(`chap5_check.py`, Aerosonde at `Va = 25`, `gamma = 0`):

| Quantity | Agreement |
|---|---|
| `a_theta1`, `a_theta2`, `a_theta3`, `a_V3` | exact to 6–7 digits |
| `a_V1`, `dT_dVa` | ~1e-4 |
| `A_lon` | 3.8e-05 scaled |
| `B_lon` | 1.4e-10 scaled |

**4. Closed loop.** A 10 m altitude step settles to under 0.5 m with airspeed
held, pitch command inside its limit and throttle inside `[0,1]`, both airframes.

### Sample output, Aerosonde at Va = 25 m/s

```
alpha* = 2.848 deg    delta_e* = -7.100 deg    delta_t* = 0.7734
residual |f(x*,u*)| = 1.9e-15

                 wn rad/s       zeta   period s
short exact       11.0168     0.4443      0.637
short approx      11.0059     0.4449      0.637
phugoid exact      0.4988     0.2858     13.144
phugoid approx     0.4752     0.3262     13.986
```

The short-period approximation is nearly exact; the phugoid approximation is
cruder, which is the expected and instructive result of section 5.6.

## Two things about the published reference

The book's companion code (`mavsim_public`) is worth being precise about, because
it is **an unfilled student skeleton, not an implementation** — every physics body
is blank, and files such as `compute_trim.m` contain `state0 =` with no
right-hand side, which will not even parse. Its README confirms that solutions
are available to instructors only. Nothing here is ported from it; the physics
comes from the text. Two findings followed from using its numbers as an oracle.

**1. The published trim point is not a force equilibrium.** Evaluating the
dynamics at its `trim_state`/`trim_input` gives `u_dot = -0.85 m/s²`. The cause:
the propeller model produces about **0.95 N** of thrust at `delta_t = 0.6768,
Va = 25`, while the drag is **10.3 N**. The true equilibrium throttle is
**0.7734**, which is what `compute_lon_trim` returns. This is consistent with
that repository's own README admitting parameter trouble and recommending
`Va = 35` for the Aerosonde.

So `check_lon` does **not** compare its trim against the published one. It
linearises at *their* point to compare `A_lon`/`B_lon` apples-to-apples, and
separately requires *its own* trim to be a genuine equilibrium. Both checks are
then tight, instead of one being confounded by the other.

**2. The Aerosonde parameters here are the Python values, not the MATLAB ones.**
That repository ships two disagreeing copies of the same aircraft:

| | its MATLAB file | its Python file (used here) |
|---|---|---|
| `C_D_0` | 0.043 | **0.0424** |
| `C_D_alpha` | 0.030 | **0.132** |
| `C_D_p` | 0.0 | 0.043 |
| `c` | 0.19 | **0.18994** |

The reference numbers were generated with the Python set. `C_D_alpha` differs by
4.4×, so this is not rounding. Verified: the Python values reproduce
`a_V1 = 0.2888454` to seven digits and the MATLAB values do not. There is a
warning comment in `lon_params.m` next to them — **do not "correct" them back.**

Also established along the way: drag is `C_D_0 + C_D_alpha*alpha`, *not* the
induced-drag polar, so `C_D_p`, aspect ratio and Oswald efficiency are genuinely
unused. And the legacy propeller model cannot be what generated the reference —
it gives `dT_dVa = -6.43` against the required `-2.35` — so only the newer
`C_T`/`C_Q` model is implemented.

## Deliberate choices

- **Finite differences default to central with a scaled step**, which is O(h²).
  A `'forward'`/`0.01` mode exists only to reproduce the published matrices,
  which were generated that way: `compute_lon_ss(x,u,P,'forward',0.01)`. That
  step is far too large for accuracy, so it is not the default. Chasing the
  reference's own discretisation error would be the wrong goal.
- **Pitch damping uses the chord**, `c/(2Va)`. Using the span instead is a real
  error that appears in some fixed-wing Simulink models and over-weights pitch
  damping by the ratio `b/c`, roughly 5× for the canard.
- **No quaternions.** With `phi = psi = v = p = r = 0` the only attitude variable
  is `theta`, an Euler angle with no singularity anywhere in the longitudinal
  envelope. The quaternion-to-Euler conversion layer that the book's companion
  code needs for its 13-state propagation is unnecessary here and is omitted.
- **The canard has no stall model.** Its `alpha0`, `M` and `epsilon` were never
  identified, so `P.stall = false` and lift stays linear. Do not encode "no
  stall" as `alpha0 = inf`; the sigmoid returns `NaN`.
- **Autopilot integrators live in a state struct, never in `persistent`.**
  `check_lon` runs the control loop more than once, and `persistent` state would
  make the second run silently wrong.
- **`theta_max` limits the pitch command, not the pitch state.** A damped
  second-order loop overshoots its command, so `|theta|` can briefly exceed the
  limit. `check_lon` asserts tightly on the command and uses a 25% band on the
  state as a divergence guard.
- **`prop_quad` thrust is not clamped at zero.** The `C_T` polynomial
  extrapolates negative at high advance ratio, which is unphysical but keeps the
  trim residual smooth and sign-changing so `fzero` stays well posed. The
  returned `delta_t` is range-checked instead.

## Not implemented

Lateral dynamics entirely (`A_lat`, `B_lat`, roll/course/sideslip loops),
Dryden wind and gusts, actuator lag, sensors and state estimation (chapters 7–8),
and the legacy propeller model. Control-surface *saturation* is kept, because
`delta_e_max` feeds the `kp_theta` formula.

## Reference

R. W. Beard and T. W. McLain, *Small Unmanned Aircraft: Theory and Practice*,
Princeton University Press, 2012. Chapters 3–6 and Appendix F.
