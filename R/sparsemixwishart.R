# WRAPPER FUNCTION
# Function for fitting Sparse Mixture of matrix Wishart for Model-based Clustering  ---------------------
#' @param data
#' @param K
#' @param penalty
#' @param penalize_diag
#' @param control
#' @param verbose
#'
#' @export
sparsemixwishart <- function(data,
                             K = 2,
                             penalty = 10,
                             penalize_diag=FALSE,
                             control = EM_controls(),
                             verbose = interactive()) {

  # The best model is the one that maximizes the BIC

  call <- match.call()
  data_dim <- dim(data)

  p <- data_dim[1]
  N <- data_dim[3]

  all_hyperparameters <-
    expand.grid(K = K,
                penalty = penalty)

  n_different_models <- nrow(all_hyperparameters)

  models_container <-
    vector(mode = "list", length = n_different_models)

  if (verbose) {
    cat("Fitting:\n")
    utils::flush.console()
    pbar <- utils::txtProgressBar(min = 0,
                                  max = n_different_models,
                                  style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }

  for (model in 1:n_different_models) {
    models_container[[model]] <-
      tryCatch(em_sparse_mix_wishart(
        data = data,
        K = all_hyperparameters[model, "K"],
        penalty = all_hyperparameters[model, "penalty"],
        penalize_diag = penalize_diag,
        control = control,
        data_dim=data_dim,
      ),error = function(e) {
        list(bic = NA)
      })

    if (verbose) {
      ipbar <- ipbar + 1
      utils::setTxtProgressBar(pbar, ipbar)
    }
  }

  models_BICS <- sapply(models_container, "[[", "bic")
  max_bic_model <- which.max(models_BICS)

  selected_model <- models_container[[max_bic_model]]
  selected_model$BIC <- cbind(all_hyperparameters, bic=models_BICS)

  if (verbose) {
    cat(
      "\nModel with K=",
      selected_model$K,
      ", penalty=",
      all_hyperparameters[max_bic_model, "penalty"],
      " returns the highest BIC.",
      sep = ""
    )
  }
  selected_model

}
