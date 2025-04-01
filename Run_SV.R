# This file assumes that Run_sampling.R has been run for the given number of rho values and repetitions
library(data.table)
source("./Estimated_SV.R")

# The parameters of the experiment
m = 10
n_train = 1000
n_test = 1000
rhos = c(0, 0.2, 0.5, 0.9)
rho_equi = TRUE
betas = c(2, 10, 0.25, -3, -1, 1.5, -0.5, 10, 1.25, 1.5, -2)
n_combinations_array = seq(4, 2^m-2, 2)
verbose_now = FALSE

# Where rds objects are saved
folder_save = "./rds"

# Create list of all feature combinations
all_coalitions = unlist(lapply(0:m, utils::combn, x = m, simplify = FALSE), recursive = FALSE)
dt_all_coalitions = data.table(features = sapply(all_coalitions, function(x) paste(x, collapse = ",")))[, id := .I]

# The number of repetitions. Change to seq(500) for 500 repetitions
repetitions = seq(10)

# Loop over the rhos and the repetitions
for (rho_idx in seq_along(rhos)) {
  # Get the current rho
  rho = rhos[rho_idx]

  # Get the true value
  file_name = paste("M", m, "n_train", n_train, "n_test", n_test,  "rho", rho, "equi", rho_equi,
                    "betas", paste(as.character(betas), collapse = "_"), sep = "_")


  # Get the model and the true Shapley values
  save_file_name_setup = file.path(folder_save, paste0(file_name, "_model.rds"))
  save_file_name_true = file.path(folder_save, paste0(file_name, "_true.rds"))

  # Load the true explanations
  message("Loading true Shapley values")
  if (!file.exists(save_file_name_true)) stop("The file with the true Shapley values do not exist.")
  true_explanations = readRDS(save_file_name_true)
  dt_vS = true_explanations$internal$output$dt_vS
  shap_names <- true_explanations$internal$parameters$feature_names
  # To support that the true Shapley values are estimated using either the old or new version of shapr
  if (is.null(true_explanations$shapley_values)) {
    dt_true_mat = as.matrix(true_explanations$shapley_values_est[,-1])
  } else {
    dt_true_mat = as.matrix(true_explanations$shapley_values[,-1])
  }
  message("Done loading true Shapley values")

  # Iterate over the repetitions
  for (repetition_idx in seq_len(repetitions)) {

    # Get the current repetition
    repetition = repetitions[repetition_idx]

    # Small printout to the user
    message(sprintf("Working on rho = %g (%d of %d) and repetition = %d (%d of %d).",
                    rho, rho_idx, length(rhos), repetition, repetition_idx, length(repetitions)))


    message("Loading presampled coalitions unique")
    presampled_coalitions_unique = file.path(folder_save, paste0("Unique_M_", m, "_repetition_", repetition, ".rds"))
    if (!file.exists(presampled_coalitions_unique)) stop("The unique sampling file does not exist for this repetition.")
    presampled_coalitions_unique = readRDS(presampled_coalitions_unique)
    message("Done loading presampled coalitions unique")

    message("Loading presampled coalitions paired")
    presampled_coalitions_paired = file.path(folder_save, paste0("Paired_M_", m, "_repetition_", repetition, ".rds"))
    if (!file.exists(presampled_coalitions_paired)) stop("The paired sampling file does not exist for this repetition.")
    presampled_coalitions_paired = readRDS(presampled_coalitions_paired)
    message("Done loading presampled coalitions paired")

    message("Loading presampled coalitions PySHAP")
    presampled_coalitions_PySHAP = file.path(folder_save, paste0("PySHAP_M_", m, "_repetition_", repetition, ".rds"))
    if (!file.exists(presampled_coalitions_PySHAP)) stop("The PySHAP sampling file does not exist for this repetition.")
    presampled_coalitions_PySHAP = readRDS(presampled_coalitions_PySHAP)
    message("Done loading presampled coalitions PySHAP")

    message("Loading presampled coalitions PySHAP*")
    presampled_coalitions_PySHAPs = file.path(folder_save, paste0("PySHAPstar_M_", m, "_repetition_", repetition, ".rds"))
    if (!file.exists(presampled_coalitions_PySHAPs)) stop("The PySHAP sampling file does not exist for this repetition.")
    presampled_coalitions_PySHAPs = readRDS(presampled_coalitions_PySHAPs)
    message("Done loading presampled coalitions PySHAP*")


    # Data.table to store the results for this repetition
    MAE_dt = data.table("Rho" = rho,
                        "Repetition" = repetition,
                        "N_S" = n_combinations_array,
                        "Unique" = NaN,
                        "Paired" = NaN,
                        "Paired C-Kernel" = NaN,
                        "PySHAP" = NaN,
                        "PySHAP*" = NaN,
                        "PySHAP* C-Kernel" = NaN)

    n_combination_idx = 50
    for (n_combination_idx in seq_along(n_combinations_array)) {

      # Get the current number of coalitions
      n_combination = n_combinations_array[n_combination_idx]

      # Small printout to the user
      message(sprintf("Working on rho = %g (%d of %d), repetition = %d (%d of %d), and n_combination = %d (%d of %d).",
                      rho, rho_idx, length(rhos), repetition, repetition_idx, length(repetitions), n_combination, n_combination_idx, length(n_combinations_array)))

      ## Unique ----------------------------------------------------------------------------------------------------------
      if (verbose_now) message("Working on Unique")
      # Get the n_combinations coalitions to include
      presampled_coalitions =
        presampled_coalitions_unique$all_coalitions[seq(presampled_coalitions_unique$dt_N_S_and_L[N_S == n_combination, L])]

      # Get the X data.table
      X_now = create_X_dt_unique_and_paired(m = m, presampled_coalitions = presampled_coalitions, dt_all_coalitions = dt_all_coalitions)

      # Compute the approximated Shapley values
      dt_kshap_unique =
        compute_SV_values(X_now = X_now, dt_all_coalitions = dt_all_coalitions, dt_vS = dt_vS, shap_names = shap_names)

      # Get the MAE between the approximated and full Shapley values
      MAE_dt[n_combination_idx, "Unique" := mean(abs(dt_true_mat - as.matrix(dt_kshap_unique[,-1])))]

      ## Paired ----------------------------------------------------------------------------------------------------------
      if (verbose_now) message("Working on Paired")
      # Get the n_combinations coalitions to include
      presampled_coalitions =
        presampled_coalitions_paired$all_coalitions[seq(presampled_coalitions_paired$dt_N_S_and_L[N_S == n_combination, L])]

      # Get the X data.table
      X_now = create_X_dt_unique_and_paired(m = m, presampled_coalitions = presampled_coalitions, dt_all_coalitions = dt_all_coalitions)

      # Compute the approximated Shapley values
      dt_kshap_paired =
        compute_SV_values(X_now = X_now, dt_all_coalitions = dt_all_coalitions, dt_vS = dt_vS, shap_names = shap_names)

      # Get the MAE between the approximated and full Shapley values
      MAE_dt[n_combination_idx, "Paired" := mean(abs(dt_true_mat - as.matrix(dt_kshap_paired[,-1])))]

      ## Paired C-kernel ----------------------------------------------------------------------------------------------------------
      if (verbose_now) message("Working on Paired C-kernel")
      # Get the n_combinations coalitions to include
      presampled_coalitions =
        presampled_coalitions_paired$all_coalitions[seq(presampled_coalitions_paired$dt_N_S_and_L[N_S == n_combination, L])]

      # Get the X data.table
      X_now = create_X_dt_unique_and_paired(m = m, presampled_coalitions = presampled_coalitions, dt_all_coalitions = dt_all_coalitions)

      # Insert the corrected Shapley kernel weights
      shapley_reweighting(X = X_now, reweight = "on_all_cond")

      # Compute the approximated Shapley values
      dt_kshap_paired_c_kernel =
        compute_SV_values(X_now = X_now, dt_all_coalitions = dt_all_coalitions, dt_vS = dt_vS, shap_names = shap_names)

      # Get the MAE between the approximated and full Shapley values
      MAE_dt[n_combination_idx, "Paired C-Kernel" := mean(abs(dt_true_mat - as.matrix(dt_kshap_paired_c_kernel[,-1])))]


      ## PySHAP ------------------------------------------------------------------------------------------------------
      if (verbose_now) message("Working on PySHAP")

      # Figure out which list to look at
      dt_id = presampled_coalitions_PySHAP$look_up$dt_n_comb_needed_sample[N_S == n_combination, dt_id]

      # Get the n_combinations coalitions to include
      to_this_index = presampled_coalitions_PySHAP$samples[[dt_id]]$dt_N_S_and_L_small[N_S == n_combination, L]
      presampled_coalitions = copy(presampled_coalitions_PySHAP$samples[[dt_id]]$all_coalitions_small[seq(to_this_index)])
      prefixed_coalitions = copy(presampled_coalitions_PySHAP$samples[[dt_id]]$dt_res)

      # Get the X data.table
      X_now = create_X_dt_PySHAP(m = m,
                                 presampled_coalitions = presampled_coalitions,
                                 prefixed_coalitions = copy(prefixed_coalitions),
                                 dt_all_coalitions = dt_all_coalitions,
                                 version_scaled = TRUE)

      # Compute the approximated Shapley values
      dt_PySHAP =
        compute_SV_values(X_now = X_now, dt_all_coalitions = dt_all_coalitions, dt_vS = dt_vS, shap_names = shap_names)

      # Get the MAE between the approximated and full Shapley values
      MAE_dt[n_combination_idx, "PySHAP" := mean(abs(dt_true_mat - as.matrix(dt_PySHAP[,-1])))]

      ## PySHAP* ------------------------------------------------------------------------------------------------------
      if (verbose_now) message("Working on PySHAP*")

      # Figure out which list to look at
      dt_id = presampled_coalitions_PySHAPs$look_up$dt_n_comb_needed_sample[N_S == n_combination, dt_id]

      # Get the n_combinations coalitions to include
      to_this_index = presampled_coalitions_PySHAPs$samples[[dt_id]]$dt_N_S_and_L_small[N_S == n_combination, L]
      presampled_coalitions = copy(presampled_coalitions_PySHAPs$samples[[dt_id]]$all_coalitions_small[seq(to_this_index)])
      prefixed_coalitions = copy(presampled_coalitions_PySHAPs$samples[[dt_id]]$dt_res)

      # Get the X data.table
      X_now = create_X_dt_PySHAP(m = m,
                                 presampled_coalitions = presampled_coalitions,
                                 prefixed_coalitions = copy(prefixed_coalitions),
                                 dt_all_coalitions = dt_all_coalitions,
                                 version_scaled = TRUE)

      # Compute the approximated Shapley values
      dt_PySHAPstar =
        compute_SV_values(X_now = X_now, dt_all_coalitions = dt_all_coalitions, dt_vS = dt_vS, shap_names = shap_names)

      # Get the MAE between the approximated and full Shapley values
      MAE_dt[n_combination_idx, "PySHAP*" := mean(abs(dt_true_mat - as.matrix(dt_PySHAPstar[,-1])))]

      ## PySHAP* C-kernel  ------------------------------------------------------------------------------------------------------
      if (verbose_now) message("Working on PySHAP* C-Kernel")

      # Figure out which list to look at
      dt_id = presampled_coalitions_PySHAPs$look_up$dt_n_comb_needed_sample[N_S == n_combination, dt_id]

      # Get the n_combinations coalitions to include
      to_this_index = presampled_coalitions_PySHAPs$samples[[dt_id]]$dt_N_S_and_L_small[N_S == n_combination, L]
      presampled_coalitions = copy(presampled_coalitions_PySHAPs$samples[[dt_id]]$all_coalitions_small[seq(to_this_index)])
      prefixed_coalitions = copy(presampled_coalitions_PySHAPs$samples[[dt_id]]$dt_res)

      # Get the X data.table
      X_now = create_X_dt_PySHAPstar(m = m,
                                     presampled_coalitions = presampled_coalitions,
                                     prefixed_coalitions = copy(prefixed_coalitions),
                                     dt_all_coalitions = dt_all_coalitions)

      # Compute the approximated Shapley values
      dt_PySHAPstar_ckernel =
        compute_SV_values(X_now = X_now, dt_all_coalitions = dt_all_coalitions, dt_vS = dt_vS, shap_names = shap_names)

      # Get the MAE between the approximated and full Shapley values
      MAE_dt[n_combination_idx, "PySHAP* C-Kernel" := mean(abs(dt_true_mat - as.matrix(dt_PySHAPstar_ckernel[,-1])))]

      if (verbose_now) print(MAE_dt[n_combination_idx,])

    } # End combinations
    print(MAE_dt)

    # Melt the data.table
    MAE_dt_long = melt(MAE_dt, id.vars = c("Rho", "Repetition", "N_S"), variable.name = "Strategy", value.name = "MAE")
    library(ggplot2)
    ggplot(MAE_dt_long, aes(x = N_S, y = MAE, col = Strategy)) +
      geom_line(linewidth = 0.65) +
      scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      )

    saveRDS(MAE_dt_long, file.path(folder_save, paste0(file_name, "_MAE_repetition_", repetition, ".rds")))

  } # End repetition
} # End rho

