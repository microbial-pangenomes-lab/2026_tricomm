# Running the Chi.bio model scripts under Octave

These MATLAB scripts were written for MATLAB, but their **numerical/model parts
run under GNU Octave** (tested with Octave 8.4) using the two small shim
functions in this folder. Only the *plotting* (`tiledlayout`, `nexttile`,
`heatmap`, `sgtitle`, `parula`) is not available in Octave — those calls error
out **after** the model has finished computing, so the numbers still come out.

## What's in this folder

| File | Purpose |
|------|---------|
| `readmatrix.m`  | Shim: Octave has no `readmatrix`; wraps `dlmread` to read the `data/*.csv` files. |
| `nlparci.m`     | Shim: reimplements MATLAB's `nlparci` (95% CIs from the Jacobian) with no Statistics toolbox. |
| `dYdt_model.m`  | Copy of the 10-state simulation ODE (`1_simulation/linear_scale_version/dYdt.m`), renamed so it never shadows the repo's own `dYdt.m` files. Used by `run_model.m`. |
| `run_model.m`   | **Parameterized model runner** — pick a parameter set, run the model, save the full time course to CSV/TSV. Plotting-free. |
| `README.md`     | This file. |

The shims work by being on Octave's *load path*. You put them there with the
`-p` (a.k.a. `--path`) flag: `octave-cli -p <folder> ...`. Nothing in the
original scripts is modified.

---

## 1. Run the model with your own parameters and save to CSV/TSV  ← main use case

`run_model.m` is a plotting-free version of
`1_simulation/linear_scale_version/tricomm_lin.m`. It integrates the 10-state
model and writes one row per time point.

**Run it from *inside* this folder** (the default data path is `../data`):

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
| `DATA_DIR` | `'../data'` | Path to the `data/` folder holding the input CSVs. |
| `OUTFILE` | `'model_output.csv'` | Output file; `.csv` → comma, `.tsv` → tab. |

**Model parameters** (the biology — the "set of parameters" to explore)

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `epsilon` | `1` | Antibiotic-induced killing rate (acts on sensitive strains). |
| `eta` | `0.005` | Conjugation rate — how fast the plasmid transfers, converting sensitive → resistant. |
| `Ab_init` | `0.18` | Initial antibiotic concentration (decays over time). |
| `scale_factor` | `1.1` | Growth-rate multiplier for resistant strains vs. sensitive (`1` = no cost/benefit; encodes plasmid fitness cost). Sets `mu_Ar = mu_A*scale_factor`, etc. |
| `mu_A`, `mu_B`, `mu_C` | `3e-4`, `5e-4`, `4e-4` | Growth rates of the sensitive A/B/C strains (cross-feeding: growth needs the two amino acids each strain consumes). |
| `alpha_A`, `alpha_B`, `alpha_C` | `0.7645`, `0.0567`, `0.1202` | Amino-acid production rates by A/B/C (A→P, B→L, C→M). |
| `delta` | `0.0692` | Common death / dilution rate for all strains. |

See the `tricomm_lin.m` header and `dYdt_model.m` for the exact equations.

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
- Uses `../../data/pOXA48_nonmLAlt_interpolated.csv` (relative path — run from
  the script's own directory as shown).
- The fit stops at `OPTIONS(14)=1000` function evaluations (set inside the
  script). Increase that value for tighter convergence — same knob as in MATLAB.

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

The data-regeneration scripts under `data/` (`*_interpolated.m`) use
`readtable` with string filtering and the raw `plate_counts_*.tsv` files; those
are **not** ported (Octave's `readtable` is limited) and aren't needed — the
interpolated CSVs they produce are already in `data/`.

---

## Requirements
- GNU Octave 8.x (`octave-cli`). No extra Octave packages required — the shims
  rely only on core functions (`ode45`, `dlmread`, `betaincinv`).
- A display is only needed for plotting; the model/CSV export runs headless.
