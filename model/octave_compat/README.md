# Running the Chi.bio model scripts under Octave

These MATLAB scripts were written for MATLAB, but their **numerical/model parts
run under GNU Octave** (tested with Octave 8.4) using the two small shim
functions in this folder. Only the *plotting* (`tiledlayout`, `nexttile`,
`heatmap`, `sgtitle`, `parula`) is not available in Octave — those calls error
out **after** the model has finished computing, so the numbers still come out.

> **Start at [`../README.md`](../README.md) instead** if you just want to run the
> model. `run_sim.m` and `run_sweep.m` take a parameter *file* as a command-line
> argument and need no editing; this file documents the shims, the published
> MATLAB scripts, and the confidence-interval fit.

## What's in this folder

| File | Purpose |
|------|---------|
| `readmatrix.m`  | Shim: Octave has no `readmatrix`; wraps `dlmread` to read the `data/model/*.csv` files. |
| `nlparci.m`     | Shim: reimplements MATLAB's `nlparci` (95% CIs from the Jacobian) with no Statistics toolbox. |
| `dYdt_model.m`  | Copy of the 10-state simulation ODE (`1_simulation/linear_scale_version/dYdt.m`), renamed so it never shadows the repo's own `dYdt.m` files. |
| `run_sim.m`     | **CLI**: run one parameter file, report diagnostics, write the time course. See [`../README.md`](../README.md) §2. |
| `run_sweep.m`   | **CLI**: sweep parameters over an N-D grid; summary table and/or all trajectories as one tidy TSV. See [`../README.md`](../README.md) §3, and §3a for the `3_grid_map` replacement. |
| `dYdt_free.m`   | The ODE the CLIs integrate: the published `dYdt` with the per-community masking dropped, so the drug, the plasmid and every compartment are set by the caller. See [`../README.md`](../README.md) §2a. |
| `tricomm_check_presets.m` | Integrates all five communities × both antibiotics through `dYdt_model.m` (the published RHS) and `dYdt_free.m` and compares — currently bit-identical. Run it after touching either file. |
| `tricomm_load_params.m` | Reads a `.par` file into a struct. Defines the file format. |
| `tricomm_simulate.m`    | Wraps the ODE: initial conditions, community overrides, `ode45`/`lsode`. |
| `tricomm_R.m`           | Coexistence ratio `R`, the `alpha_C` root for `R = 1`, equilibrium composition. |
| `tricomm_mse.m`         | Fit of any 7-parameter set to the KAN no-drug replicates. |
| `tricomm_fate.m`, `tricomm_grid.m`, `tricomm_apply_sets.m`, `tricomm_argparse.m`, `tricomm_opt.m` | Small helpers for the two CLIs. |
| `run_model.m`   | The earlier runner — parameters live in a CONFIG block you edit or override inline. Still works; `run_sim.m` supersedes it. |
| `README.md`     | This file. |

The shims work by being on Octave's *load path*. You put them there with the
`-p` (a.k.a. `--path`) flag: `octave-cli -p <folder> ...`. Nothing in the
original scripts is modified.

---

## 1. Run the model with your own parameters and save to CSV/TSV

> Superseded by `run_sim.m` — see [`../README.md`](../README.md) §2, which does
> everything below plus parameter files, diagnostics and sweeps. `run_model.m` is
> kept because existing notes and scripts refer to it, and it still works.

`run_model.m` is a plotting-free version of
`1_simulation/linear_scale_version/tricomm_lin.m`. It integrates the 10-state
model and writes one row per time point.

**Run it from *inside* this folder** (the default data path is `../../data/model`):

```bash
cd octave_compat
octave-cli -p . --eval "run_model"
```

That writes `model_output.csv` with these columns:

```
time, As, Ar, Bs, Br, Cs, Cr, P, L, M, Ab, A_tot, B_tot, C_tot, A_CFUmL, B_CFUmL, C_CFUmL
```
- `As..Ab` — the 10 model states (A/B/C sensitive & resistant, amino acids P/L/M, antibiotic Ab)
- `*_tot`   — sensitive + resistant, in model units
- `*_CFUmL` — the totals × 1e7 (CFU/mL, as plotted in the paper).  A=LM, B=PM, C=PL.

### Choosing the parameter set — two ways

**(a) Edit the CONFIG block** at the top of `run_model.m` (antibiotic, community,
`epsilon`, `eta`, `Ab_init`, `scale_factor`, the `mu_*` / `alpha_*` / `delta`
rates, `TMAX`), then re-run `octave-cli -p . --eval "run_model"`.

**(b) Override inline** without editing the file — set any variables before
calling `run_model` in the same `--eval`:

```bash
# no-drug community, custom transfer rate, 1000 h, tab-separated output
octave-cli -p . --eval "community=0; eta=0.02; TMAX=1000; OUTFILE='nodrug.tsv'; run_model"

# sweep a parameter from the shell
for e in 0.2 0.5 1.0; do
  octave-cli -p . --eval "epsilon=$e; OUTFILE='eps_$e.csv'; run_model"
done
```

Output format follows the `OUTFILE` extension: **`.csv` → comma, `.tsv` → tab**.

### Parameters you can change

All of these are in the CONFIG block of `run_model.m` (or can be set inline).
Strains: **A = LM, B = PM, C = PL**. Amino acids in the model: **P, L, M**.

**Experiment / run setup**

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `antibiotic` | `1` | `1` = ampicillin 16×MIC, `2` = colistin (then set `dose`). Picks which data columns are used and the antibiotic decay rate `gamma`. |
| `dose` | `2` | Only for `antibiotic=2`: `1` = 1×MIC, `2` = 2×MIC. |
| `community` | `3` | `0` = no-drug, `1..4` = drug communities. **`0` forces `epsilon=eta=Ab_init=0`; `1` and `4` force `eta=0`** (as in the original model). Also sets the initial conditions. |
| `TMAX` | `[]` | Final time in hours. `[]` = last experimental time point (306 h here). |
| `DATA_DIR` | `'../../data/model'` | Path to the `data/model/` folder holding the input CSVs. |
| `OUTFILE` | `'model_output.csv'` | Output file; `.csv` → comma, `.tsv` → tab. |

**Model parameters** (the biology — the "set of parameters" to explore)

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `epsilon` | `1` | Antibiotic-induced killing rate (acts on sensitive strains). |
| `eta` | `0.005` | Conjugation rate — how fast the plasmid transfers, converting sensitive → resistant. |
| `Ab_init` | `0.18` | Initial antibiotic concentration (decays over time). |
| `scale_factor` | `1.1` | Growth-rate multiplier for resistant strains vs. sensitive (`1` = no cost/benefit; encodes plasmid fitness cost). Sets `mu_Ar = mu_A*scale_factor`, etc. |
| `mu_A`, `mu_B`, `mu_C` | `0.0003`, `0.0005`, `0.0004` | Second-order growth-rate **constants** (units 1/(h·conc²) — it is `mu_A*L*M` that has units of 1/h, so `mu_A` on its own is not a growth rate). Fitted; see section 2. |
| `alpha_A`, `alpha_B`, `alpha_C` | `0.7645`, `0.0567`, `0.1202` | Amino-acid production rates by A/B/C (A→P, B→L, C→M). Fitted; see section 2. |
| `delta` | `0.0692` | Common death / dilution rate for all strains. Fitted. Also sets *when* the community first dips (dip time ≈ const/`delta`). |

See the `tricomm_lin.m` header and `dYdt_model.m` for the exact equations.

### Model form and where the parameters come from

The equations are the **original published mass-action model** — `dYdt.m` is
unchanged in every copy — and the parameters are the **original published set**
(`ori_par`), fitted by the authors against the pOXA48 no-drug mean. Nothing in
this repo's model code deviates from what was published.

Alternatives were explored and are documented, not wired in:

- two collaborator refits against the KAN replicates (`calc_par` with µ locked to
  the measured isolated ratio, and `fit_par` with µ free);
- a set that keeps `calc_par` but **solves** `alpha_C` from the coexistence
  condition, which is the only one of the four that does not grow without bound.

All four are ready to run as parameter files in `model/params/`; see
**`model/README.md`** for the runner and **`model/PARAMETER_SET_ASSESSMENT.md`**
for the comparison that led to keeping `ori_par`.

---

## 2. Parameter estimation + 95% confidence intervals

`6_heatmap_alpha_ic/ci_calculation/confidence_interval.m` runs the SIMPLEXL
(Nelder–Mead) fit, then computes CIs. It has **no plotting**, so it runs to
completion and prints a table.

```bash
cd 6_heatmap_alpha_ic/ci_calculation
octave-cli -p ../../octave_compat --eval "confidence_interval"
```

Notes:
- Fits the **pOXA48 no-drug mean** — columns 2–4 of
  `data/model/pOXA48_nonmLAlt_interpolated.csv`, 27 interpolated points. That file
  is built from an *older* snapshot of the plate counts
  (`data/model/plate_counts_pOXA48_old.tsv`), not from the current
  `data/pOXA48/`; see `model/README.md` §6. Keep that snapshot, or the fit's input
  is gone.
- The committed data path is `../../data/pOXA48_nonmLAlt_interpolated.csv`, which
  does not resolve in this repo layout. Point it at
  `../../../data/model/` before running.
- **Seven free parameters** (all three `mu`, all three `alpha`, `delta`) and a
  **linear** squared-error objective on the three strains — not log10. Because the
  strains span ~3 orders of magnitude, this objective is dominated by the most
  abundant one.
- **Running it as committed does not return `ori_par`.** It stops on its iteration
  cap (`options(14) = 1000`; Octave prints a "maximum iterations reached" warning)
  at `mu ≈ (4.9e-4, 1.1e-3, 5.6e-4)`, `alpha ≈ (0.584, 0.0520, 0.0888)`,
  `delta ≈ 0.0585`.
- **`ori_par` is not a converged optimum of this objective either.** Its MSE is
  14.24, and restarting Nelder–Mead from it walks steadily down to 10.51 at
  `alpha_A ≈ 0.54`, `delta ≈ 0.0571`. `ori_par` is a snapshot of a partially
  converged run, and the published CIs are the linearisation *at that snapshot* —
  their midpoints are exactly `ori_par`. So the CIs reproduce; the optimisation
  that produced the point estimate does not.
- Recomputing those CIs at `ori_par` under Octave gives the same centres but
  intervals **1.5–1.8× wider** than published (e.g. `alpha_A` 0.4597–1.0693
  against 0.5932–0.9359). The widening is not a constant factor, so it is not a
  degrees-of-freedom convention; it points to a slightly different residual
  vector, most likely MATLAB's `ode45` against Octave's.
- The published `alpha` CIs are recorded in `6_heatmap_alpha_ic/heatmap_alpha_ic.m`
  (cases 1 and 3), which is what that script sweeps. `mu` is held at `ori_par` in
  all three cases.
- An interior steady state requires
  `alpha_A*alpha_B*alpha_C/delta^3 == (alpha_A+alpha_B+alpha_C)/delta + 2` — the
  growth rates cancel out entirely. With `a = alpha_A/delta` etc., the ratio
  `R = abc/(a+b+c+2)` is > 1 for a growing community, < 1 for a decaying one, and
  = 1 for a neutral line of equilibria; its distance from 1 sets the rate.
  `ori_par` sits at `R = 1.0077`, so it grows slowly and without bound. Report it
  for any set with `run_sim.m` (see `model/README.md` §5).
- `OPTIONS(14)` (inside the script) caps the Nelder–Mead iterations — same knob
  as in MATLAB; increase for tighter convergence.

---

## 3. The original simulation / sweep / heatmap scripts

These all **compute correctly** but **fail at the plotting stage** in Octave.
Run them from their own directory with this folder on the path:

```bash
cd 1_simulation/linear_scale_version
octave-cli -p ../../octave_compat --eval "tricomm_lin"   # integrates, then errors at the first plot call
```

The same is true for `3_grid_map/grid_map.m`, `4_heatmap_BrCr/heatmap_BrCr.m`,
`5_heatmap_alpha_1/heatmap_alpha_1.m`, `6_heatmap_alpha_ic/heatmap_alpha_ic.m`:
the `ode45` sweeps finish and the result matrices (e.g. `crossing_time_1`) are
in the workspace before the error.

- **If you only need the numbers**, prefer `run_model.m` (section 1) or add your
  own `csvwrite`/`dlmwrite` before the plotting section.
- **If you need the figures**, run those scripts in real MATLAB or
  MATLAB Online (https://matlab.mathworks.com) — Octave does not implement
  `tiledlayout`/`heatmap`.

The data-regeneration scripts under `data/model/` (`*_interpolated.m`) use
`readtable` with string filtering and the raw `plate_counts.tsv` files (in the
per-experiment folders `data/pOXA48/`, `data/PN23/`; the `_old` variant lives in
`data/model/`); those are **not** ported (Octave's `readtable` is limited) and
aren't needed — the interpolated CSVs they produce are already in `data/model/`.

---

## Requirements
- GNU Octave 8.x (`octave-cli`). No extra Octave packages required — the shims
  rely only on core functions (`ode45`, `dlmread`, `betaincinv`).
- A display is only needed for plotting; the model/CSV export runs headless.
