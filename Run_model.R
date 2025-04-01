# File to create the models to explain
library(mvtnorm)
library(xgboost)
library(tidymodels)
library(data.table)
library(future)
library(progressr)
library(mgcv)

# Where rds objects are saved
folder_save = "./rds"

# The number of features
m = 10

# The number of training and test samples/explicands
n_train = 1000
n_test = 1000

# The correlation levels
rhos = c(0.0, 0.2, 0.5, 0.9)
rho_equi = TRUE

# The coefficients in the data generating process
betas = c(2, 10, 0.25, -3, -1, 1.5, -0.5, 10, 1.25, 1.5, -2)

# Mean of the multivariate Gaussian distribution
mu = rep(0, times = m)

# We use the Gaussian approach
approach = "gaussian"

# The number of Monte Carlo samples used to estimate the contribution functions.
# Larger values means that we obtain Shapley values closer to the true values,
# but also drastically increases the running time
n_MC_samples = 10000

# Number of cores to run the Shapley values computations on
n_workers = 4

# Iterate over the rhos
rho_idx = 1
for (rho_idx in seq_along(rhos)) {
  # Get the current rho
  rho = rhos[rho_idx]

  # Small message to the user
  message(paste0("Working on rho = ", rho, " (", rho_idx, " of ", length(rhos), ")"))

  # Create the covariance matrix
  sigma = matrix(ncol = m, nrow = m)
  if (rho_equi) {
    sigma = matrix(rho, ncol = m, nrow = m)
  } else {
    for (i in seq(1, m-1)) {
      for (j in seq(i+1, m))
        sigma[i,j] = sigma[j,i] = rho^abs(i-j)
    }
  }
  diag(sigma) = 1

  # Make some of the save file names
  file_name = paste("M", m, "n_train", n_train, "n_test", n_test,  "rho", rho, "equi", rho_equi,
                    "betas", paste(as.character(betas), collapse = "_"), sep = "_")
  save_file_name_setup = file.path(folder_save, paste0(file_name, "_model.rds"))
  save_file_name_true = file.path(folder_save, paste0(file_name, "_true.rds"))

  # Create the data and predictive model
  message("Create the data and predictive model.")
  gamma_coefficients = c(1, -2, 2, -3) # Pair-wise interactions between [X1, X2], [X3,X5], [X4,X8]
  gamma_terms = list(c("X1", "X2"), c("X3", "X5"), c("X4", "X8"), c("X9", "X10"))
  eta_coefficients = c(3, -1, -2) # [X1, X3, X7] and [X2, X6, X8] and [X3, X8, X10]
  eta_terms = list(c("X1", "X3", "X7"), c("X2", "X6", "X8"), c("X3", "X8", "X10"))
  alpha_coefficients = c(4) # [X1, X4, X7, X9]
  alpha_terms = list(c("X1", "X4", "X7", "X9"))

  # Set seed for reproducibility
  seed_setup = 2000
  set.seed(seed_setup)

  # Make Gaussian data
  data_train = data.table(rmvnorm(n = n_train + 10, mean = mu, sigma = sigma))
  data_test  = data.table(rmvnorm(n = n_test + 10,  mean = mu, sigma = sigma))
  colnames(data_train) = paste("X", seq(m), sep = "")
  colnames(data_test) = paste("X", seq(m), sep = "")

  # Make the response
  response_train_org = as.vector(cbind(1, as.matrix(data_train)) %*% betas)
  response_test_org = as.vector(cbind(1, as.matrix(data_test)) %*% betas)

  tmp_add_interaction_values = function(terms, coefficients, data) {
    rowSums(sapply(seq_along(terms), function(i) {
      terms_now = terms[[i]]
      coefficient_now = coefficients[i]
      return(apply(as.matrix(data[,..terms_now]), 1, prod) * coefficient_now)
    }))
  }

  response_train = response_train_org +
    tmp_add_interaction_values(gamma_terms, gamma_coefficients, data_train) +
    tmp_add_interaction_values(eta_terms, eta_coefficients, data_train) +
    tmp_add_interaction_values(alpha_terms, alpha_coefficients, data_train)

  response_test = response_test_org +
    tmp_add_interaction_values(gamma_terms, gamma_coefficients, data_test) +
    tmp_add_interaction_values(eta_terms, eta_coefficients, data_test) +
    tmp_add_interaction_values(alpha_terms, alpha_coefficients, data_test)

  # Remove the ten largest absolute values as this DGP can create some outliers
  train_indices = order(abs(response_train), decreasing = TRUE)[-c(1:10)]
  test_indices = order(abs(response_test), decreasing = TRUE)[-c(1:10)]
  response_train = response_train[train_indices]
  response_train_org = response_train_org[train_indices]
  response_test = response_test[test_indices]
  response_test_org = response_test_org[test_indices]
  data_test = data_test[test_indices]
  data_train = data_train[train_indices]

  # # Plot the data
  # plot_data = data.table(values = c(response_train_org, response_train), type = rep(c("without interactions", "with interactions"), times = c(n_train, n_train)))
  # ggplot(plot_data, aes(values, fill = type)) + geom_histogram(alpha = 0.5, aes(y = after_stat(density)), position = 'identity')
  #
  # plot_data = data.table(values = c(response_test_org, response_test), type = rep(c("without interactions", "with interactions"), times = c(n_test, n_test)))
  # ggplot(plot_data, aes(values, fill = type)) + geom_histogram(alpha = 0.5, aes(y = after_stat(density)), position = 'identity')

  # Put together the data
  data_train_with_response = copy(data_train)[,y := response_train]
  data_test_with_response  = copy(data_test)[,y := response_test]

  # Predictive model: use tidymodels to do CV to find best hyperparameters
  regression.workflow = workflows::add_recipe(
    workflows::add_model(
      workflows::workflow(),
      parsnip::boost_tree(trees = hardhat:::tune(),
                          tree_depth = hardhat:::tune(),
                          learn_rate = hardhat:::tune(),
                          engine = "xgboost",
                          mode = "regression")),
    recipes::recipe(as.formula("y ~ ."),
                    data = data_train_with_response))
  regression.results <- tune::tune_grid(
    object = regression.workflow,
    resamples = rsample::vfold_cv(data = data_train_with_response, v = 5),
    grid = expand.grid(tree_depth = c(2, 4, 6, 8),
                       trees = c(5, 10, 25, 50, 100, 200, 250, 500, 1000, 1500, 2000),
                       learn_rate = c(0.05, 0.1, 0.2)),
    metrics = yardstick::metric_set(yardstick::rmse),
    control = tune::control_grid(verbose = TRUE)
  )
  print(tune::show_best(regression.results, metric = "rmse", n = 10))

  # Look at the accuracy of the model compared to a basic lm and a GAM model
  predictive_model = lm(y ~ ., data = data_train_with_response)
  message(sprintf("LM: Training MSE = %g and test MSE = %g.",
                  mean((predict(predictive_model, data_train_with_response) - data_train_with_response$y)^2),
                  mean((predict(predictive_model, data_test_with_response) - data_test_with_response$y)^2)))

  predictive_model = gam(as.formula(paste0("y ~ ", paste0("ti(X", seq(m), ")", collapse = " + "))),
                         data = data_train_with_response)
  message(sprintf("GAM: Training MSE = %g and test MSE = %g.",
                  mean((predict(predictive_model, data_train_with_response) - data_train_with_response$y)^2),
                  mean((predict(predictive_model, data_test_with_response) - data_test_with_response$y)^2)))

  # Train an xgboost model with the best hyperparameters
  best_results = tune::select_best(regression.results, metric = "rmse")
  predictive_model = xgboost(data = as.matrix(data_train), label = response_train,
                             nrounds = best_results$trees,
                             params = list("eta" = best_results$learn_rate, "max_depth" = best_results$tree_depth),
                             verbose = FALSE)

  message(sprintf("XGBOOST: Training MSE = %g and test MSE = %g.",
                  mean((predict(predictive_model, as.matrix(data_train)) - data_train_with_response$y)^2),
                  mean((predict(predictive_model, as.matrix(data_test)) - data_test_with_response$y)^2)))

  # plot(predict(predictive_model, as.matrix(data_train)), data_train_with_response$y)
  # plot(predict(predictive_model, as.matrix(data_test)), data_test_with_response$y)

  # Get the prediction zero, i.e., the phi0 Shapley value.
  prediction_zero = mean(response_train)

  # Collect the variables
  save_list = list(seed_setup = seed_setup,
                   betas = betas,
                   data_train = data_train,
                   data_test = data_test,
                   response_train = response_train,
                   response_test = response_test,
                   data_train_with_response = data_train_with_response,
                   data_test_with_response = data_test_with_response,
                   prediction_zero = prediction_zero,
                   predictive_model = predictive_model)

  # Save the results
  message("Start saving the model.")
  if (file.exists(save_file_name_setup)) stop(paste0("The file `", save_file_name_setup, "` already exists."))
  saveRDS(save_list, save_file_name_setup)
  message("Saved the model.")

  # Compute the true explanations
  message(paste0("Start computing the true explanations (n_workers = ", n_workers, ")."))

  # Set future in the right plan
  if (n_workers > 1) {
    future::plan(multisession, workers = n_workers)
  } else {
    future::plan(sequential)
  }

  # Specify the progressr bar
  progressr::handlers("cli")

  # We are either using a non-Gaussian approach or a non-linear predictive model
  # Compute the true explanations
  progressr::with_progress({
    true_explanations <- shapr::explain(
      model = predictive_model,
      x_explain = data_test,
      x_train = data_train,
      approach = approach,
      phi0 = prediction_zero,
      exact = TRUE,
      iterative = FALSE,
      n_MC_samples = n_MC_samples,
      gaussian.mu = mu,
      gaussian.cov_mat = sigma,
      seed = 1
    )
    }, enable = TRUE)

  # Set future back to sequential plan
  future::plan(sequential)

  # Save the true explanations just in case
  message("Start saving the true explanations.")
  if (file.exists(save_file_name_true)) stop(paste0("The file `", save_file_name_true, "` already exists."))
  saveRDS(true_explanations, save_file_name_true)
  message("Saved the true explanations.")
}
