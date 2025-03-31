# Libraries -------------------------------------------------------------------------------------------------------
library(data.table)
library(stringr)
source("./shapr_functions.R")

# Functions -------------------------------------------------------------------------------------------------------
compute_SV_values <- function(X_now, dt_all_coalitions, dt_vS, shap_names) {
  # Get the weight matrix
  W_now <- weight_matrix(X = X_now, normalize_W_weights = TRUE)

  # Use the pre-computed v(S) data and only extract the relevant rows (combinations)
  dt_vS_now <- as.matrix(dt_vS[X_now[, id_combination_full], -"id_combination"])

  # Compute the new Shapley values
  dt_kshap <- data.table::as.data.table(t(W_now %*% dt_vS_now))
  colnames(dt_kshap) <- c("none", shap_names)
  return(dt_kshap)
}


create_X_dt_unique_and_paired <- function(m, presampled_coalitions, dt_all_coalitions, weight_zero_m = 10^6) {
  # Find weights for given number of features
  n_features <- seq(m - 1)
  n <- sapply(n_features, choose, n = m)
  w <- shapley_weights(m = m, N = n, n_features) * n
  p <- w / sum(w)

  # String version
  # Insert all sampled coalitions into a data table and find their frequencies
  dt_freq <- data.table::data.table(features = presampled_coalitions)[, .(shapley_weight = .N), by = features]

  # Get the number of features in each coalition
  dt_freq[, n_features := stringr::str_count(features, ",") + 1]

  # Add the number of coalitions of each size
  dt_freq[, N := n[n_features]]
  dt_freq[, p := p[n_features]]

  # Get the id_combination if we had used all combinations
  dt_freq[, id_combination_full := dt_all_coalitions[dt_freq, id, on = "features"]]

  # Convert from string to list of integer vectors. stringr is faster than base and stringi
  dt_freq[, features := lapply(stringr::str_split(features, ","), as.integer)]

  # Add the empty and grand coalitions
  dt_freq <- rbindlist(
    list(
      data.table(
        features = list(integer(0)),
        shapley_weight = weight_zero_m,
        n_features = 0L,
        N = 1L,
        p = NA,
        id_combination_full = 1
      ),
      dt_freq,
      data.table(
        features = list(1:m),
        shapley_weight = weight_zero_m,
        n_features = as.integer(m),
        N = 1L,
        p = NA,
        id_combination_full = 2^m
      )
    )
  )
  data.table::setorder(dt_freq, "id_combination_full")
  dt_freq[, id_combination := .I]
  data.table::setcolorder(
    dt_freq,
    c("id_combination", "id_combination_full", "features", "n_features", "N", "shapley_weight", "p")
  )

  # Optional to match the old setup
  dt_freq[, N := as.integer(N)]
  dt_freq[, shapley_weight := as.integer(shapley_weight)]
  dt_freq[, n_features := as.integer(n_features)]

  return(dt_freq)
}

create_X_dt_PySHAP <- function(m,
                               presampled_coalitions,
                               prefixed_coalitions,
                               dt_all_coalitions,
                               weight_zero_m = 10^6,
                               version_scaled = TRUE) {

  # Find weights for given number of features
  n_features <- seq(m - 1)
  n <- sapply(n_features, choose, n = m)
  w <- shapley_weights(m = m, N = n, n_features) * n
  p <- w / sum(w)

  # Weight the different coalition sizes (PySHAP version)
  num_subset_sizes <- as.integer(ceiling((m - 1) / 2))
  num_paired_subset_sizes <- as.integer(floor((m - 1) / 2))
  weight_vector <- sapply(seq(num_subset_sizes), function(i) (m - 1.0) / (i * (m - i)))
  weight_vector[seq(num_paired_subset_sizes)] <- 2 * weight_vector[seq(num_paired_subset_sizes)]
  weight_vector <- weight_vector / sum(weight_vector)

  # Get the sampled coalitions for which we have to compute the frequencies
  if (!is.null(prefixed_coalitions)) {
    presampled_coal_wo_prefixed_coal <- presampled_coalitions[-seq(nrow(prefixed_coalitions))]
  } else {
    presampled_coal_wo_prefixed_coal <- presampled_coalitions
  }

  # String version
  # Insert all sampled coalitions into a data table and find their frequencies
  dt_freq <- data.table::data.table(features = presampled_coal_wo_prefixed_coal)[, .(shapley_weight = .N), by = features]

  # Fix the weights according to the technique in PySHAP
  if (version_scaled) {
    if (is.null(prefixed_coalitions)) {
      num_full_subsets <- 0
      weight_left <- sum(weight_vector)
    } else {
      num_full_subsets <- length(prefixed_coalitions[.N - 1, features][[1]]) # This relies on the list version
      weight_left <- sum(weight_vector[-seq(num_full_subsets)])
    }
    dt_freq[, shapley_weight := shapley_weight * weight_left / sum(shapley_weight)]
  }

  # Convert the list column to a comma-separated string for each row
  if (!is.null(prefixed_coalitions)) {
    prefixed_coalitions[, features := sapply(features, function(x) paste(unlist(x), collapse = ","))]
    setnames(prefixed_coalitions, "w", "shapley_weight")
  }

  # Put together with the prefixed samples
  dt_freq <- rbind(prefixed_coalitions, dt_freq)

  # Get the number of features in each coalition
  dt_freq[, n_features := stringr::str_count(features, ",") + 1]

  # Add the number of coalitions of each size
  dt_freq[, N := n[n_features]]
  dt_freq[, p := p[n_features]]

  # Get the id_combination if we had used all combinations
  dt_freq[, id_combination_full := dt_all_coalitions[dt_freq, id, on = "features"]]

  # Convert from string to list of integer vectors. stringr is faster than base and stringi
  dt_freq[, features := lapply(stringr::str_split(features, ","), as.integer)]

  # Add the empty and grand coalitions
  dt_freq <- rbindlist(
    list(
      data.table(features = list(integer(0)), shapley_weight = weight_zero_m, n_features = 0L, N = 1L, p = NA, id_combination_full = 1),
      dt_freq,
      data.table(features = list(1:m), shapley_weight = weight_zero_m, n_features = as.integer(m), N = 1L, p = NA, id_combination_full = 2^m)
    )
  )
  data.table::setorder(dt_freq, "id_combination_full")
  dt_freq[, id_combination := .I]
  data.table::setcolorder(dt_freq, c("id_combination", "id_combination_full", "features", "n_features", "N", "shapley_weight", "p"))
  # dt_freq

  # Optional to match the old setup
  dt_freq[, N := as.integer(N)]
  dt_freq[, n_features := as.integer(n_features)]

  # plot(dt_freq[-c(1, .N), id_combination_full], dt_freq[-c(1, .N), shapley_weight])
  #
  # dt_freq[, sum(shapley_weight), by = N][-1]
  return(dt_freq)
}

create_X_dt_PySHAPstar <- function(m,
                                         presampled_coalitions,
                                         prefixed_coalitions,
                                         dt_all_coalitions,
                                         weight_zero_m = 10^6) {
  # Find weights for given number of features
  n_features <- seq(m - 1)
  n <- sapply(n_features, choose, n = m)
  w <- shapley_weights(m = m, N = n, n_features) * n
  p <- w / sum(w)

  # Weight the different coalition sizes (PySHAP version)
  num_subset_sizes <- as.integer(ceiling((m - 1) / 2))
  num_paired_subset_sizes <- as.integer(floor((m - 1) / 2))
  weight_vector <- sapply(seq(num_subset_sizes), function(i) (m - 1.0) / (i * (m - i)))
  weight_vector[seq(num_paired_subset_sizes)] <- 2 * weight_vector[seq(num_paired_subset_sizes)]
  weight_vector <- weight_vector / sum(weight_vector)

  # Get the sampled coalitions for which we have to compute the frequencies
  if (!is.null(prefixed_coalitions)) {
    presampled_coal_wo_prefixed_coal <- presampled_coalitions[-seq(nrow(prefixed_coalitions))]
    num_full_subsets <- length(prefixed_coalitions[.N - 1, features][[1]]) # This relies on the list version
    weights_remaining <- weight_vector[-seq(num_full_subsets)]
    weight_left <- sum(weights_remaining)
  } else {
    presampled_coal_wo_prefixed_coal <- presampled_coalitions
    num_full_subsets <- 0
    weights_remaining <- weight_vector
    weight_left <- sum(weight_vector)
  }

  # Get the number of coalitions of all sizes
  n_coal_of_each_size <- choose(m, seq(m - 1))

  # Get the number of coalitions in the sizes that are not deterministically included
  if (num_full_subsets >= floor(m / 2)) stop("Too many full subsets. No sampling is done.")
  n_coal_of_each_size_reamaining <- n_coal_of_each_size[seq(num_full_subsets + 1, m - 1 - num_full_subsets)]

  # Get the shapley kernel weights for the remaining coalition sizes
  p_reamaining <- p[seq(num_full_subsets + 1, m - 1 - num_full_subsets)]
  p_reamaining <- p_reamaining / sum(p_reamaining)

  # Get the shapley kernel weight for each coalition
  shapley_kernel_weight_reweighted <- p_reamaining / n_coal_of_each_size_reamaining

  # Pad it such that index corresponds to caolition size
  shapley_kernel_weight_reweighted <- c(rep(0, num_full_subsets), shapley_kernel_weight_reweighted, rep(0, num_full_subsets))

  # Convert the list column to a comma-separated string for each row
  if (!is.null(prefixed_coalitions)) {
    prefixed_coalitions[, features := sapply(features, function(x) paste(unlist(x), collapse = ","))]
    setnames(prefixed_coalitions, "w", "shapley_weight")

    # Add that these coalitions are NOT sampled
    prefixed_coalitions[, sampled := FALSE]
  }

  # String version
  # Insert all sampled coalitions into a data table and find their frequencies
  dt_freq <- data.table::data.table(features = presampled_coal_wo_prefixed_coal)[, .(shapley_weight = .N), by = features]

  # Add that these coalitions are sampled
  dt_freq[, sampled := TRUE]

  # Put together with the prefixed samples
  dt_freq <- rbind(prefixed_coalitions, dt_freq)

  # Get the number of features in each coalition
  dt_freq[, n_features := stringr::str_count(features, ",") + 1]

  # Add the number of coalitions of each size
  dt_freq[, N := n[n_features]]
  dt_freq[, p := p[n_features]]

  # Get the id_combination if we had used all combinations and sort it
  dt_freq[, id_combination_full := dt_all_coalitions[dt_freq, id, on = "features"]]
  data.table::setorder(dt_freq, "id_combination_full")

  # Get the number of samples it took the sample the `dt_freq[sampled == TRUE, .N]` unique
  # coalitions for those coalitions that are not deterministically included.
  K <- dt_freq[sampled == TRUE, sum(shapley_weight)]

  # Ensure that the shapley weight column is numeric, as it is int when prefixed_coalitions is NULL
  dt_freq[, shapley_weight := as.numeric(shapley_weight)]

  # Add the reweighted Shapley kernel weights
  dt_freq[sampled == TRUE, shapley_weight := shapley_kernel_weight_reweighted[n_features]]

  # Compute the corrected values shapley weights (see paper)
  dt_freq[sampled == TRUE, shapley_weight := 2 * shapley_weight / (1 - (1 - 2 * shapley_weight)^(K / 2))]

  # Reweight the weights such that they sums to the remaining weight
  dt_freq[sampled == TRUE, shapley_weight := weight_left * shapley_weight / sum(shapley_weight)]

  # Convert from string to list of integer vectors. stringr is faster than base and stringi
  dt_freq[, features := lapply(stringr::str_split(features, ","), as.integer)]

  # Remove the sampled column as we no longer need it
  dt_freq[, sampled := NULL]

  # Add the empty and grand coalitions
  dt_freq <- rbindlist(
    list(
      data.table(features = list(integer(0)), shapley_weight = weight_zero_m, n_features = 0L, N = 1L, p = NA, id_combination_full = 1),
      dt_freq,
      data.table(features = list(1:m), shapley_weight = weight_zero_m, n_features = as.integer(m), N = 1L, p = NA, id_combination_full = 2^m)
    )
  )

  # Create new id column and reorder the columns
  dt_freq[, id_combination := .I]
  data.table::setcolorder(dt_freq, c("id_combination", "id_combination_full", "features", "n_features", "N", "shapley_weight", "p"))

  # Optional to match the old setup
  dt_freq[, N := as.integer(N)]
  dt_freq[, n_features := as.integer(n_features)]

  return(dt_freq)
}
