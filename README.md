# simple-sequential-testing-experiments
Simulation experiments for the paper


- `run_simulations.R` - the main R script that runs all of the simulation experiments
- `methods` - a folder with the implementation of the sequential tests under comparison

To reproduce the results simply run `run_simulations.R`. It will print tables with FDR (false detection rate), power, and sample savings for experiments with both normal and non-normal data.