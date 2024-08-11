# YEt Another Sequential Test
Experiments for the paper


- `run_simulations.R` - an R script that runs the simulation experiments
- `evaluate_on_online_retail.R` - an R script that validates the method using a public real-world dataset (Online Retail)
- `evaluate_on_real_world_data_2.R` - an R script that validates the method using real-world data from a major e-commerce platform 
- `data_generation.R` - an R script with classes for data generation
- `methods` - a folder with the implementation of the sequential tests under comparison
- `utils.R` - an R script with helper classes and functions

To reproduce the simulation results run `run_simulations.R`. It will print tables with FDR (false detection rate), power, and sample savings for experiments with both normal and non-normal data.

To reproduce the results of the validation on the [Online Retail](https://archive.ics.uci.edu/dataset/352/online+retail) dataset, run `evaluate_on_online_retail.R`.
