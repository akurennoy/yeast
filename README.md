# simple-sequential-testing-experiments
Simulation experiments for the paper


- `20230321_sim_blogg.R` - main experiment with normal increments
- `sim_t.R` and `sim_gamma.R` - experiments with increments simulated using t and Gamma distributions, respectively

To reproduce the result of an experiment simply run the respective R script. It will print tables with FDR (false detection rate), power, and sample savings. The number of replications is set to 100000, so script execution will take some time.

The results presented in the main body of the paper can be found in the "stream" column of the respective tables. It corresponds to the continuous monitoring regime.
