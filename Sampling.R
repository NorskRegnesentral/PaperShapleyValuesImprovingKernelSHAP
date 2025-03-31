# Libraries -------------------------------------------------------------------------------------------------------
library(data.table)
source("./shapr_functions.R")

# Functions -------------------------------------------------------------------------------------------------------
coalition_sampling_unique = function(m, n_combinations = 2^m - 2,  n_sample_scale = 5, return_coalitions = FALSE,
                                     seed = NULL, verbose = TRUE) {
  if (n_combinations > 2^m - 2) stop("n_combinations is larger than 2^m.")
  if (!is.null(seed)) set.seed(seed)

  # Find weights for given number of features
  n_features <- seq(m - 1)
  n <- sapply(n_features, choose, n = m)
  w <- shapley_weights(m = m, N = n, n_features) * n
  p <- w / sum(w)

  # List to store all the sampled coalitions
  all_coalitions = c()

  # Variable to keep track of the number of unique coalitions
  unique_coalitions = 0

  # Variable to keep track of the iteration number
  iteration = 1

  # Loop until we have enough unique samples
  while (unique_coalitions < n_combinations) {

    # Sample the coalition sizes
    n_features_sample <- sample(
      x = n_features,
      size = n_sample_scale*n_combinations,
      replace = TRUE,
      prob = p
    )

    # Sample the coalitions
    coalitions <- sample_features_cpp(m, n_features_sample)

    # Convert the coalitions to strings such that we can compare them
    coalitions = sapply(coalitions, paste, collapse = ",")

    # Add the new coalitions to the previously sampled coalitions
    all_coalitions = c(all_coalitions, coalitions)

    # Get the cumulative number of unique coalitions for each coalition in all_coalitions
    dt_cumsum = data.table(coalitions = all_coalitions, N_S = cumsum(!duplicated(all_coalitions)))[, L := .I]

    # Extract rows where the N_S value increases (i.e., where we sample a new unique coalition)
    dt_N_S_and_L <- dt_cumsum[N_S != shift(N_S, type = "lag", fill = 0)]

    # Get the number of unique coalitions
    unique_coalitions = dt_N_S_and_L[.N, N_S]

    # Message to user
    if (verbose) {
      message(paste0("Iteration ", iteration, ": N_S = ", unique_coalitions,
                     ", Sampled = ", n_sample_scale*n_combinations*iteration, "."))
    }

    # Update the iteration number
    iteration = iteration + 1
  }

  # Post processing: keep only the coalitions until n_combinations
  all_coalitions = all_coalitions[seq(dt_N_S_and_L[N_S == n_combinations, L])]
  if (length(unique(all_coalitions)) != n_combinations) stop("Not the right number of unique coalitions")

  # Return
  if (return_coalitions) {
    return(list(dt_N_S_and_L = dt_N_S_and_L, all_coalitions = all_coalitions))
  } else {
    return(dt_N_S_and_L)
  }
}

coalition_sampling_paired = function(m,
                                     n_combinations = 2^m - 2,
                                     n_sample_scale = 5,
                                     return_coalitions = FALSE,
                                     seed = NULL,
                                     verbose = TRUE) {

  if (n_combinations > 2^m - 2) stop("n_combinations is larger than 2^m.")
  if (!is.null(seed)) set.seed(seed)

  # Find weights for given number of features
  n_features <- seq(m - 1)
  n <- sapply(n_features, choose, n = m)
  w <- shapley_weights(m = m, N = n, n_features) * n
  p <- w / sum(w)

  # List to store all the sampled coalitions
  all_coalitions = c()

  # Variable to keep track of the number of unique coalitions
  unique_coalitions = 0

  # Variable to keep track of the iteration number
  iteration = 1

  # Loop until we have enough unique samples
  while (unique_coalitions < n_combinations) {

    # Sample the coalition sizes
    message("Getting the coalition sizes")
    n_features_sample <- sample(
      x = n_features,
      size = n_sample_scale*n_combinations,
      replace = TRUE,
      prob = p
    )

    # Sample the coalitions
    message("Getting coalitions")
    feature_sample <- sample_features_cpp(m, n_features_sample)

    # Get the paired coalitions
    message("Making the paired")
    feature_sample_paired <- lapply(feature_sample, function(x, m) {seq(m)[-x]}, m = m)

    # Merge the coalitions in alternating fashion as we do paired sampling (i.e., first is S and second is Sbar and so on)
    coalitions = c(rbind(feature_sample, feature_sample_paired))

    # Convert the coalitions to strings such that we can compare them
    message("Converting to strings")
    coalitions = sapply(coalitions, paste, collapse = ",")

    # Add the new coalitions to the previously sampled coalitions
    all_coalitions = c(all_coalitions, coalitions)

    # Get the cumulative number of unique coalitions for each coalition in all_coalitions
    message("Getting cumsum")
    dt_cumsum = data.table(coalitions = all_coalitions, N_S = cumsum(!duplicated(all_coalitions)))[, L := .I]

    # Extract rows where the N_S value increases (i.e., where we sample a new unique coalition)
    message("Getting shift")
    dt_N_S_and_L <- dt_cumsum[N_S != shift(N_S, type = "lag", fill = 0)]

    # Get the number of unique coalitions
    unique_coalitions = dt_N_S_and_L[.N, N_S]

    # Message to user
    if (verbose) {
      message(paste0("Iteration ", iteration, ": N_S = ", unique_coalitions,
                     ", Sampled = ", n_sample_scale*n_combinations*iteration, "."))
    }

    # Update the iteration number
    iteration = iteration + 1
  }

  # Post processing: keep only the coalitions until n_combinations
  all_coalitions = all_coalitions[seq(dt_N_S_and_L[N_S == n_combinations, L])]
  if (length(unique(all_coalitions)) != n_combinations) stop("Not the right number of unique coalitions")

  # Return
  if (return_coalitions) {
    return(list(dt_N_S_and_L = dt_N_S_and_L, all_coalitions = all_coalitions))
  } else {
    return(dt_N_S_and_L)
  }
}



coalition_sampling_kernelSHAP = function(m,
                                         n_combinations = 2^m - 2,
                                         n_sample_scale = 3,
                                         return_coalitions = TRUE,
                                         seed = NULL,
                                         verbose = TRUE,
                                         always_pair_coalitions = TRUE) {
  if (n_combinations > 2^m - 2) stop("n_combinations is larger than 2^m.")
  if (!is.null(seed)) set.seed(seed)

  # Number of features
  M = m

  # weight the different coalition sizes
  num_subset_sizes = as.integer(ceiling((M - 1) / 2))
  num_paired_subset_sizes = as.integer(floor((M - 1) / 2))
  weight_vector = sapply(seq(num_subset_sizes), function(i) (M - 1.0) / (i * (M - i)))
  weight_vector[seq(num_paired_subset_sizes)] = 2*weight_vector[seq(num_paired_subset_sizes)]
  weight_vector = weight_vector / sum(weight_vector)

  # Variable that will store the normalized probability of sampling the remaining colaition sizes
  remaining_weight_vector = copy(weight_vector)

  # Array to store the number of combinations needed to include the different coalition sizes
  n_comb_needed = NULL

  # Find the number of combinations needed to include the different coalition sizes
  subset_size = 1
  for (subset_size in seq(num_subset_sizes)) {

    # Get the number of (paired) coalitions of this subset size
    nsubsets = choose(M, subset_size)
    if (subset_size <= num_paired_subset_sizes) nsubsets = 2 * nsubsets

    # Get the expected number of samples needed to sample nsubsets coalitions of size
    # `subset_size` using the normalized sampling probability
    n_comb_needed_now = ceiling(nsubsets / remaining_weight_vector[subset_size])

    # Add the number of coalitions of smaller sizes that are included
    if (subset_size > 1) n_comb_needed_now = n_comb_needed_now + 2 * sum(choose(M, seq(subset_size - 1)))

    # Store the new values
    n_comb_needed = c(n_comb_needed, n_comb_needed_now)

    # Update the probability of the remaining coalition sizes such that they sum to 1
    if (remaining_weight_vector[subset_size] < 1.0) {
      remaining_weight_vector = remaining_weight_vector / (1 - remaining_weight_vector[subset_size])
    }
  }

  # Create a data table with max number of coalitions before we include the smaller coalition size. Ensure even numbers
  n_comb_needed = sapply(n_comb_needed, function(x) ifelse(x %% 2 == 0, x - 2, x - 1))
  n_comb_needed[n_comb_needed >= n_combinations] = n_combinations
  n_comb_needed[length(n_comb_needed)] = n_combinations
  dt_n_comb_needed = data.table(dt_id = seq_along(n_comb_needed), N_S = n_comb_needed)
  dt_n_comb_needed[, N_S_fixed := fifelse(dt_id == 1, 0, 2 * sapply(dt_id, function(id) sum(choose(M, seq_len(id - 1)))))]
  dt_n_comb_needed

  # Create a look up table.
  # The idea now is that if for each value of N_S, we can get which result list to look at by looking at
  # `dt_n_comb_needed_sample[N_S == 916, dt_id]`.
  dt_n_comb_needed_sample = data.table(N_S = seq(2, n_combinations, 2),
                                       dt_id = sapply(seq(2, n_combinations, 2), function(x) which.max(n_comb_needed >= x)))


  id_now_idx = 2
  id_max = length(dt_n_comb_needed$dt_id)
  full_res = lapply(seq_along(dt_n_comb_needed$dt_id), function(id_now_idx) {
    id_now = dt_n_comb_needed$dt_id[id_now_idx]

    # Get the number of unique coalitions to sample
    N_S_now = dt_n_comb_needed[dt_id == id_now, N_S]

    # data table to store the coalitions that are pre-defined to be included with the corresponding normalized Shapley kernel weights
    dt_res = NULL

    # For all id_now larger than 1, we include all coalitions of sizes less than `id_now`
    if (id_now > 1) {
      subset_size = 1
      for (subset_size in seq(id_now - 1)) {
        feature_sample = unlist(lapply(subset_size, utils::combn, x = M, simplify = FALSE), recursive = FALSE)
        w = weight_vector[subset_size] / choose(M, subset_size)
        if (subset_size <= num_paired_subset_sizes) {
          # Add paired sampled and half the weight
          feature_sample = c(rbind(feature_sample, lapply(feature_sample, function(x, M) {seq(M)[-x]}, M = M)))
          w = w / 2
        }
        dt_res_now = data.table(features = feature_sample, w = w)
        dt_res = rbind(dt_res, dt_res_now)
      }
    }

    # Convert to string
    dt_res_features_string = sapply(dt_res$features, paste, collapse = ",")

    # Then we need to sample the remaining features
    # add random samples from what is left of the subset space
    nfixed_samples = ifelse(is.null(dt_res), 0, nrow(dt_res))
    samples_left = N_S_now - nfixed_samples

    num_full_subsets = id_now - 1

    if (num_full_subsets != num_subset_sizes) {
      # Get the normalized weights for the remaining coalition sizes
      remaining_weight_vector = copy(weight_vector)
      if (always_pair_coalitions) {
        remaining_weight_vector = remaining_weight_vector / 2 # because we draw two samples each below
      } else {
        remaining_weight_vector[seq(num_paired_subset_sizes)] = remaining_weight_vector[seq(num_paired_subset_sizes)] / 2 # because we draw two samples each below
      }
      if (num_full_subsets > 0) remaining_weight_vector = remaining_weight_vector[-seq(num_full_subsets)] # Remove the fully sampled coalition size
      remaining_weight_vector = remaining_weight_vector / sum(remaining_weight_vector)


      # List to store all the sampled coalitions
      all_coalitions = c()

      # Variable to keep track of the number of unique coalitions
      unique_coalitions = nfixed_samples

      # Variable to keep track of the iteration number
      iteration = 1

      # Loop until we have enough unique samples
      while (unique_coalitions < N_S_now) {

        # Sample the coalition sizes
        message(paste0("(", id_now, "/", id_max, ") ", "Getting the coalition sizes"))
        n_features_sample <- sample(
          x = length(remaining_weight_vector),
          size = n_sample_scale * samples_left,
          replace = TRUE,
          prob = remaining_weight_vector
        ) + num_full_subsets # Add the num_full_subsets to get the correct coal sizes

        # Sample the coalitions
        message(paste0("(", id_now, "/", id_max, ") ", "Getting coalitions"))
        feature_sample <- sample_features_cpp(m, n_features_sample)

        # Get the paired coalitions
        message(paste0("(", id_now, "/", id_max, ") ", "Making the paired"))
        if (always_pair_coalitions) {
          # Get the paired coalitions
          feature_sample_paired <- lapply(feature_sample, function(x, m) {seq(m)[-x]}, m = m)

          # Merge the coalitions in alternating fashion as we do paired sampling (i.e., first is S and second is Sbar and so on)
          coalitions = c(rbind(feature_sample, feature_sample_paired))

        } else {
          # In python SHAP, they do not pair the coalition of M/2 for M even.
          # This is strange as we then no longer can garante that both S and Sbar are sampled
          coalitions <- unlist(lapply(feature_sample, function(x, m) {
            if (length(x) == M / 2) {
              return(list(x))
            } else {
              return(list(x, seq(m)[-x]))
            }
          }, m = m), recursive = FALSE)
        }

        # Convert the coalitions to strings such that we can compare them
        message(paste0("(", id_now, "/", id_max, ") ", "Converting to strings"))
        coalitions = sapply(coalitions, paste, collapse = ",")

        # Add the new coalitions to the previously sampled coalitions
        all_coalitions = c(all_coalitions, coalitions)

        # Add the fixed coalitions
        if (nfixed_samples > 0) {
          all_coalitions_added = c(dt_res_features_string, all_coalitions)
        } else {
          all_coalitions_added = all_coalitions
        }

        # Get the cumulative number of unique coalitions for each coalition in all_coalitions_added
        message(paste0("(", id_now, "/", id_max, ") ", "Getting cumsum"))
        dt_cumsum = data.table(coalitions = all_coalitions_added, N_S = cumsum(!duplicated(all_coalitions_added)))[, L := .I]
        dt_cumsum

        # Extract rows where the N_S value increases (i.e., where we sample a new unique coalition)
        message(paste0("(", id_now, "/", id_max, ") ", "Getting shift"))
        dt_N_S_and_L <- dt_cumsum[N_S != shift(N_S, type = "lag", fill = 0)]

        # Get the number of unique coalitions
        unique_coalitions = dt_N_S_and_L[.N, N_S]

        if (unique_coalitions < N_S_now) stop("Not enough sampels")

        # Message to user
        if (verbose) {
          message(paste0("Iteration ", iteration, ": N_S_now = ", N_S_now,
                         ": N_S = ", unique_coalitions, ", Sampled = ", nrow(dt_cumsum), "."))
        }

        # Update the iteration number
        iteration = iteration + 1
      }

      # Stop at next limit
      dt_N_S_and_L[N_S == N_S_now]
      dt_N_S_and_L_small = dt_N_S_and_L[N_S <= N_S_now]
      all_coalitions_small = all_coalitions_added[seq(dt_N_S_and_L_small[.N, L])]
    }

    # Return
    if (return_coalitions) {
      return(list(dt_res = dt_res, dt_N_S_and_L_small = dt_N_S_and_L_small, all_coalitions_small = all_coalitions_small))
    } else {
      return(dt_N_S_and_L)
    }
  })


  final_list = list(look_up = list(dt_n_comb_needed = dt_n_comb_needed, dt_n_comb_needed_sample = dt_n_comb_needed_sample),
                    samples = full_res)

  #print(object.size(final_list), units = "MB")
  return(final_list)
}


