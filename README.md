# 2026_tricomm

Data, model and analysis code for the three-member community experiments.

**You don't need to run anything to see the results.** The processed tables
the scripts produce are committed to the repo, and every notebook is saved
with its outputs (plots included). Only re-run things if you want to reproduce
or change them.

## Layout

| Path | Contents |
|------|----------|
| `data/<experiment>/` | Raw measurements (Chi.bio logs, plate counts, flow cytometry, plate reader) for each experiment (`pOXA48`, `PN23`, `long`, `control`, `batch`, `calibration`, `pairs`), plus the notebook that parses them into tidy tables (`chibio.tsv.gz`, `plate_counts.tsv`, `flow_cytometer.tsv`, ...) |
| `data/pairs/` | Chi.bio run with every pair of strains and the full three-strain community; parsed by `analysis.ipynb` into `chibio.tsv.gz` |
| `data/model/` | Model outputs (`*.tsv`) and the interpolated experimental input to the model (`*_interpolated.csv`) |
| `model/` | MATLAB model scripts; `model/octave_compat/` has Octave-compatible versions and the `run_sim.m`/`run_sweep.m` command-line tools |
| `model/params/params.par` | Model parameters |
| `models.sh` | Runs all the simulations used in the paper and writes them to `data/model/` |
| `notebooks/` | Analysis and figure notebooks; they read from `data/` |
| `*.ipynb` (top level) | Exploratory notebooks (example runs, correlation between measurements) |
| `telegram_bot/` | Bot used to monitor the Chi.bio during experiments; not needed for the analysis |

## Environment

```bash
conda env create -f environment.yml
conda activate 2026_tricomm
```

The environment includes GNU Octave, which runs the model (tested with 8.4).

## Reproducing the outputs

Run the steps in this order; each one overwrites the committed files.

1. Parse the raw data (run from inside each experiment folder):

   ```bash
   cd data/pOXA48
   jupyter nbconvert --to notebook --execute --inplace parsing.ipynb
   ```

   Same for `data/PN23/parsing.ipynb`, `data/long/parsing.ipynb`,
   `data/control/analysis.ipynb`, `data/pairs/analysis.ipynb`,
   `data/batch/parse_time_series.ipynb` and
   `data/calibration/calibration.ipynb`.

2. Build the model input from the pOXA48 plate counts (from the repo root):

   ```bash
   python data/pOXA48/build_pOXA48_nonmLAlt_interpolated.py
   ```

3. Run the model simulations (from the repo root; this takes a while):

   ```bash
   bash models.sh
   ```

4. Re-run the figure notebooks (from inside `notebooks/`):

   ```bash
   cd notebooks
   jupyter nbconvert --to notebook --execute --inplace *.ipynb
   ```

   Figures are written as PNG/SVG next to the notebooks.
