# Source files ----------------------------------------------------------------------------------------------------
source("./Sampling.R") # Take some time as it will compile the C++ code

# Run code --------------------------------------------------------------------------------------------------------
# Note that the running time and memory consumption exponentially increases with the number of features.
# User get more messages with lower values of n_sample_scale, but then the code takes longer to execute.
# The code also run faster if n_combinations is lower than 2^m - 2, as sampling the last remanining coalitions
# take a long time, especially for the unique approach.

# The number of features and coalitions
m = 10
n_combinations = 2^m - 2

# The multiplicative number of coalitions to sample
n_sample_scale = 20

# Where to save the MAE values
folder_save = "./rds"

# The number of repetitions
repetitions = seq(100)

# The different sampling versions
versions = c("Unique", "Paired", "PySHAP", "PySHAP*")

# Loop over the repetitions and strategies
for (repetition in repetitions) {
  for (version in versions) {
    message(paste0("\nWorking on version '", version, "' and repetition ", repetition, "."))

    # Generate the coalitions
    if (version == "Unique") {
      sampled_coalitions = coalition_sampling_unique(
        m = m,
        n_combinations = n_combinations,
        n_sample_scale = n_sample_scale,
        return_coalitions = TRUE,
        seed = repetition + 1
      )
    } else if (version == "Paired") {
      sampled_coalitions = coalition_sampling_paired(
        m = m,
        n_combinations = n_combinations,
        n_sample_scale = n_sample_scale,
        return_coalitions = TRUE,
        seed = repetition + 1
      )
    } else if (version == "PySHAP") {
      sampled_coalitions = coalition_sampling_kernelSHAP(
        m = m,
        n_combinations = n_combinations,
        n_sample_scale = n_sample_scale,
        return_coalitions = TRUE,
        seed = repetition + 1,
        always_pair_coalitions = FALSE
      )
    } else if (version == "PySHAP*") {
      sampled_coalitions = coalition_sampling_kernelSHAP(
        m = m,
        n_combinations = n_combinations,
        n_sample_scale = n_sample_scale,
        return_coalitions = TRUE,
        seed = repetition + 1,
        always_pair_coalitions = TRUE
      )
      version = "PySHAPstar" # AS file on windows cannot contain *
    } else {
      stop("Version is not supported.")
    }

    # Print the size
    print(object.size(sampled_coalitions), units = "MB")

    # Save the file
    saveRDS(sampled_coalitions, file.path(folder_save, paste0(version, "_M_", m, "_repetition_", repetition, ".rds")))
  }
}
