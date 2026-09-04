# YEt Another Sequential Test

Experiments for the YEAST paper.

## Layout

| path | contents |
|---|---|
| `run_simulations.R` | the simulation experiments (normal, Student's *t*, Gamma) |
| `evaluate_on_online_retail.R` | validation on the public [Online Retail](https://archive.ics.uci.edu/dataset/352/online+retail) dataset |
| `evaluate_on_real_world_data_2.R` | validation on real-world data from an e-commerce platform |
| `data_generation.R` | data-generating classes |
| `methods/` | the sequential tests under comparison |
| `utils.R` | helper classes and functions |
| `scripts/` | environment setup, smoke test, table and figure generation |
| `tests/` | formula-level tests |

## Reproducing the results

Run everything from the repository root; the scripts `source()` each other by
relative path.

```sh
Rscript scripts/setup_environment.R   # install packages from a pinned CRAN snapshot
Rscript tests/run_all.R               # formula-level tests, seconds
Rscript scripts/smoke_test.R          # 2,000 null replications, minutes
Rscript run_simulations.R             # the full experiment, 100,000 replications
Rscript scripts/make_tables.R         # writes tables/*.tex
Rscript scripts/make_figures.R        # writes power-4.png
Rscript scripts/session_info.R > ENVIRONMENT.md
```

`run_simulations.R` writes `normal.csv`, `student.csv` and `gamma.csv`. Every
table and figure in the manuscript is generated from those files by
`scripts/make_tables.R` and `scripts/make_figures.R`; no reported figure should
be transcribed by hand. Set `YEAST_NUM_REPLICATIONS` to shorten a run while
developing.

`methods/ld_obf.R` generates `obf_bounds_100.csv` on first use, so a clean
checkout needs no pre-existing artefacts.

## Conventions

These matter when comparing the code with the manuscript.

**Sign.** The simulation code tracks treatment minus control, so a positive
relative effect places the alternative in the upper tail that the one-sided
boundary monitors. The main text of the paper uses the opposite convention, in
which the control group carries the positive sign. Appendix B of the paper is
written in the code's convention and says so.

**Aggregation unit.** The cluster-robust variance estimates use the subject
(user) as the cluster. The real-world evaluation aggregates events to the
user-month before clustering, not to the user-hour.

**Boundaries.** Every boundary-based test implements `boundary(n)`, which
returns the monitored indices and the alerting boundary at those indices on the
cumulative-sum scale; `monitor()` is inherited and is exactly
`trajectory[index] > value`. `tests/test_boundaries.R` pins each boundary to the
closed form it is meant to implement.

**Savings.** `get_savings()` reports the fraction of the planned event budget
left unused. A replication with no detection contributes exactly zero, and the
reported averages are over all replications, so they mix detection frequency
with detection earliness and are not the savings conditional on detection.
Savings are measured in events, not calendar time; the two coincide only under a
stable arrival rate.

## Benchmark implementations

The benchmarks are run on repeated-user event streams that lie outside the
settings their published guarantees cover. This is deliberate — it is the
operational question a practitioner faces — and the resulting numbers are
robustness observations about these implementations at these operating points,
not evidence about the underlying theorems. Specifically:

- `methods/msprt.R` is the Robbins normal mixture with a fixed cluster-scale
  plug-in for the increment standard deviation, at level `alpha* = 2 * alpha`.
- `methods/gavi.R` is the same boundary family: the cumulative form of the mSPRT
  boundary at `phi` is identical to the GAVI boundary at `rho = phi`. The two
  differ only in how the tuning value is specified. `tests/test_msprt_boundary.R`
  asserts this identity.
- `methods/ld_obf.R` holds the 100-look O'Brien-Fleming critical values constant
  between scheduled looks so that the boundary can be applied after every event.
  Discrete-monitoring alpha-spending baselines are compared separately, as GST
  variants, in the appendix experiment.
- `methods/sec_c_2st_ons_qda.R` applies SeqC2ST-QDA to the raw event stream.
