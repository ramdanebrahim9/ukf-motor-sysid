# ukf-motor-sysid

Adaptive Unscented Kalman Filter for DC motor state estimation — 4-state augmented filter running on real experimental data captured at 20 kHz from a custom STM32H755 acquisition setup.

---

## What this does

Estimates 4 states simultaneously from encoder ticks and current sensor readings:

| State | Symbol | Unit | Description |
|-------|--------|------|-------------|
| Position | `θ` | rad | Output shaft angle (gear ratio N=13.7335) |
| Velocity | `ω` | rad/s | Angular velocity, filtered from 64 PPR encoder |
| Current | `i` | A | Armature current, sensed via INA240 |
| Disturbance | `d` | A | Lumped model error — friction residual + unmodeled dynamics |

The filter runs **augmented**: `d` is estimated online as a 4th state with a first-order decay model (`τ_d = 0.1 s`), corrected at each step by a clipping strategy driven by the RK4 open-loop model drift.

---

## Experimental data

`EXP_172_raw.mat` contains real DC motor signals captured at **Ts = 50 µs (20 kHz)** using:

- **STM32H755ZI** dual-core (CM7 @ 480 MHz for control + UKF, CM4 for SD logging)
- **TIM1** center-aligned PWM at 20 kHz — ADC triggered at midpoint to avoid switching noise
- **INA240** current sense amplifier with differential RC anti-aliasing filter (~1.2 kHz cutoff)
- **OPA197 / OPA333** signal conditioning buffers on voltage sense lines
- **Quadrature encoder** — 64 PPR, read via hardware timer, timestamped deterministically via DWT cycle counter
- **SDMMC + FatFS** with ping-pong SRAM4 buffers for lossless high-rate logging

The experiment sweeps 12 voltage levels (1V → 12V) across ~78 seconds, each partition lasting ~6.5 s. Motor parameters were identified beforehand via a full grey-box sysid pipeline (R, L, Ke, J, B, Stribeck friction).

**Motor parameters:**

| Parameter | Value | Unit |
|-----------|-------|------|
| R | 1.438907 | Ω |
| L | 0.415 | mH |
| Ke | 0.132926 | V·s/rad |
| J | 1.263×10⁻⁴ | kg·m² |
| Gear ratio N | 13.7335 | — |
| Encoder | 64 | PPR |

---

## Filter design

### Process model — augmented ODE (`motor_ode_Aug.m`)

```
dθ/dt = ω
dω/dt = (Ke/J)·i_rk4 − friction(ω)          ← i injected from open-loop RK4
di/dt = (1/L)·(V − R·i − Ke·ω_rk4 + d)      ← ω injected from open-loop RK4
dd/dt = −d / τ_d
```

Friction model: Stribeck + viscous
```
f(ω) = tanh(ω/ε)·[tc/J + (ts/J − tc/J)·exp(−(ω/ωs)²)] + (B/J)·ω
```

RK4 injection decouples the current and mechanical states — avoids the instability of a fully coupled 4-state sigma point propagation at 20 kHz.

### Adaptive Q — `get_Q.m`

Process noise on `i` is voltage-scheduled via a lookup table (10 breakpoints, 1V→10V). `θ` and `ω` are fixed. `d` uses a constant `q_d = 0.01`.

### Adaptive R — `get_R.m`

Measurement noise is dual-scheduled:

- **`r_i`** — interpolated from a 10-point ω-indexed table (characterized from ADC noise vs operating point)
- **`r_θ`, `r_ω`** — tick-event scheduled:
  - On tick: `r_θ = 1e-6`, `r_ω = k1(V)` (trust the encoder)
  - Between ticks: `r_θ = tick_rad`, `r_ω = k2(V)` (coast with high uncertainty)

### Disturbance clipping

At each step, `d` is hard-clipped against the RK4 model drift:
```
drift = mean_I_200Hz − i_rk4
d ← clip(d, −|drift|, +|drift|)   if |drift| > thresh(V)
```
Threshold is voltage-dependent, derived from the offset characterization table.

---

## Repo structure

```
ukf-motor-sysid/
│
├── data/
│   └── EXP_172_raw.mat          # real STM32 capture — 20 kHz, ~78 s
│
├── ukf/
│   ├── compute_weights.m        # sigma point weights (alpha, beta, kappa)
│   ├── compute_sigma_points.m   # 2n+1 sigma points via Cholesky
│   ├── unscented_transform.m    # weighted mean + covariance
│   ├── ukf_predict.m            # predict step
│   └── ukf_correct.m            # correct step
│
├── model/
│   ├── rk4_step.m               # RK4 integrator — 2-state (ω, i)
│   ├── rk4_step_Aug.m           # RK4 integrator — 4-state augmented
│   ├── motor_ode_sim.m          # standard motor ODE (Stribeck friction)
│   └── motor_ode_Aug.m          # augmented ODE with disturbance state d
│
├── tuning/
│   ├── get_Q.m                  # adaptive process noise — V-scheduled
│   └── get_R.m                  # adaptive measurement noise — ω + tick scheduled
│
├── Deploy_UKF.m                 # main script — runs all regions, plots results
├── UKF_Aug_ALL_regions.pdf      # exported results — code + figures
├── .gitignore
└── README.md
```

---

## Running

```matlab
% From repo root in MATLAB:
addpath(fullfile(pwd, 'ukf'));
addpath(fullfile(pwd, 'model'));
addpath(fullfile(pwd, 'tuning'));
addpath(fullfile(pwd, 'data'));

run('Deploy_UKF.m')
```

Edit the `REGIONS` cell at the top of `Deploy_UKF.m` to select voltage partitions:

```matlab
REGIONS = {'2v', '5v', '7.5v'};   % any subset of the 12 available partitions
```

Available partitions: `1v, 2v, 3v, 4v, 5v, 6.5v, 7.5v, 8.5v, 9v, 10v, 11v, 12v`

Each region produces one figure with a 2×2 grid: **Omega · Current · Disturbance · Model Drift**

---

## Dependencies

Standard MATLAB toolboxes only — Signal Processing Toolbox (`butter`, `filtfilt`). No additional toolboxes required.

---

## Author

**Ibrahim Ramdane** — M1 VISTA (Vision, Signal, Trajectographie et Automatique), Université de Toulon  
Control systems · State estimation · Embedded real-time implementation
