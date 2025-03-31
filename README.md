# PaperShapleyValuesImprovingKernelSHAP
Relevant code to the "Improving the Weighting Strategy in KernelSHAP" paper by Lars Henry Berge Olsen and Martin Jullum.

There is a three step procedure to replicate Figure 2 in the paper.
1. Run the file `Run_sampling.R` with the desired number of repetitions `B` to generate files that generates and store the sampled coalitions using the different sampling strategies.
2. Run the file `Run_SV.R` to compute the MAE error between the estimated and true Shapley values for different number of used unique coalitions, and store the results as data tables to disk.
3. Run the file `Run_plot.R` to load the data tables with the MAE and make the plot.

Note. We have only made the models and true Shapley values available for the `M=10` simulation study, as the files are several GBs for the `M = 20` setting and it thus not suitable for GitHub.
