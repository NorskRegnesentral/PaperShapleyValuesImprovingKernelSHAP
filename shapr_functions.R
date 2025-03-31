# Libraries -------------------------------------------------------------------------------------------------------
library(Rcpp)

# Functions -------------------------------------------------------------------------------------------------------
# The functions are taken from the shapr package

Rcpp::cppFunction(
  " List sample_features_cpp(int m, IntegerVector n_features) {
    int n = n_features.length();
    List l(n);
    for (int i = 0; i < n; i++) {
        int s = n_features[i];
        IntegerVector k = sample(m, s);
        std::sort(k.begin(), k.end());
        l[i] = k;

    }
    return l;
}")

Rcpp::cppFunction("
#include <RcppArmadillo.h>
arma::mat weight_matrix_cpp(List subsets, int m, int n, NumericVector w){

    // Define objects
    int n_elements;
    IntegerVector subset_vec;
    arma::mat Z(n, m + 1, arma::fill::zeros), X(n, m + 1, arma::fill::zeros);
    arma::mat R(m + 1, n, arma::fill::zeros);

    // Populate Z
    for (int i = 0; i < n; i++) {

        // Set all elements in the first column equal to 1
        Z(i, 0) = 1;

        // Extract subsets
        subset_vec = subsets[i];
        n_elements = subset_vec.length();
        if (n_elements > 0) {
            for (int j = 0; j < n_elements; j++)
                Z(i, subset_vec[j]) = 1;
        }
    }

    // Populate X
    for (int i = 0; i < n; i++) {

        for (int j = 0; j < Z.n_cols; j++) {

            X(i, j) = w[i] * Z(i, j);
        }
    }
    R = solve(X.t() * Z, X.t());
    return R;
}", depends = c("RcppArmadillo"))

shapley_weights <- function(m, N, n_components, weight_zero_m = 10^6) {
  x <- (m - 1) / (N * n_components * (m - n_components))
  x[!is.finite(x)] <- weight_zero_m
  x
}

sum_shapley_weights <- function(m){
  coal_samp_vec <- seq(m - 1)
  n <- sapply(coal_samp_vec, choose, n = m)
  w <- shapley_weights(m = m, N = n, coal_samp_vec) * n
  return(sum(w))
}

shapley_reweighting <- function(X, reweight = "on_N", K = NULL, coal_size_fixed = NULL) {
  # Updates the Shapley weights in X based on the reweighting strategy BY REFERENCE

  if (reweight == "on_N") {
    X[-c(1,.N), shapley_weight := mean(shapley_weight), by = N]
  } else if (reweight == "on_coal_size") {
    X[-c(1,.N), shapley_weight := mean(shapley_weight), by = n_features]
  } else if (reweight == "on_all") {
    m <- X[.N, n_features]
    X[-c(1,.N), shapley_weight := shapr:::shapley_weights(m = m, N = N, n_components = n_features, weight_zero_m = 10^6)/sum_shapley_weights(m)]
  } else if (reweight == "on_N_sum") {
    X[, shapley_weight := sum(shapley_weight), by = N]
  } else if (reweight == "on_all_cond") {
    X[, shapley_weight := as.numeric(shapley_weight)]
    m <- X[.N, n_features]
    if (is.null(K)) K <- X[-c(1,.N), sum(shapley_weight)]
    X[-c(1,.N), shapley_weight := shapr:::shapley_weights(m = m, N = N, n_components = n_features, weight_zero_m = 10^6)/sum_shapley_weights(m)]
    X[-c(1,.N), cond := 1-(1-shapley_weight)^K]
    X[-c(1,.N), shapley_weight := shapley_weight/cond]
  } else if (reweight == "on_all_cond_paired") {
    X[, shapley_weight := as.numeric(shapley_weight)]
    m <- X[.N, n_features]
    if (is.null(K)) K <- X[-c(1,.N), sum(shapley_weight)]
    X[-c(1,.N), shapley_weight := shapr:::shapley_weights(m = m, N = N, n_components = n_features, weight_zero_m = 10^6)/sum_shapley_weights(m)]
    X[-c(1,.N), cond := 1-(1-2*shapley_weight)^(K/2)]
    X[-c(1,.N), shapley_weight := 2*shapley_weight/cond]
  } else if (reweight == "comb") {
    # Very ad-hoc: Half on_coal_size and half of on_all
    X[, shapley_weight1 := mean(shapley_weight), by = n_features]
    X[-c(1,.N), shapley_weight1 := shapley_weight1/sum(shapley_weight1)]
    m <- X[.N, n_features]
    X[-c(1,.N), shapley_weight2 := shapr:::shapley_weights(m = m, N = N, n_components = n_features, weight_zero_m = 10^6)/sum_shapley_weights(m)]
    X[-c(1,.N), shapley_weight2 := shapley_weight2/sum(shapley_weight2)]
    X[-c(1,.N), shapley_weight := (shapley_weight1+shapley_weight2)/2]
  }
  # strategy= "none" or something else do nothing
  return(NULL)
}

weight_matrix <- function(X, normalize_W_weights = TRUE) {
  # Fetch weights
  w <- X[["shapley_weight"]]

  if (normalize_W_weights) {
    w[-c(1, length(w))] <- w[-c(1, length(w))] / sum(w[-c(1, length(w))])
  }

  W <- weight_matrix_cpp(
    subsets = X[["features"]],
    m = X[.N][["n_features"]],
    n = X[, .N],
    w = w
  )
  return(W)
}
