#' @export
EM_controls <- function(tol = c(1e-05),
                        max_iter = 1e03,
                        type_start = c("hc", "random"),
                        linkage_hc_start= "ward.D",
                        # n_subset_start = NULL,
                        n_random_start = 50,
                        nu_max=NULL)
  # EM control parameters
{
  list(tol = tol,
       max_iter = max_iter,
       type_start = match.arg( type_start, choices = eval(formals(EM_controls)$type_start) ),
       linkage_hc_start= linkage_hc_start,
       # n_subset_start = n_subset_start,
       n_random_start = n_random_start,
       nu_max=nu_max)
}


nu_k_root_finder_equation <- function(data, p, z_k, n_k, nu_k, Sigma_inv_k) {
  # Equation (5) of \cite{Hidot2010} Pattern Recognition letter
  sum(z_k * apply(data, 3, function(s)
    determinant(x = (s %*% Sigma_inv_k) / 2, logarithm = TRUE)$modulus)) - n_k * sum(digamma(x = (nu_k -
                                                                                                    (1:p) + 1) / 2))
}

weighted_sample_S_calculator <- function(data, z, p, n_K) {
  K <- ncol(z)
  weighted_data <- vector(mode = "list", length = K)
  S_weighted <- array(dim = c(p, p, K))
  for (k in 1:K) {
    # weighted_data[[k]] <- sweep(x = data,MARGIN = 3,STATS = z[,k],FUN = "*")
    weighted_data <- sweep(
      x = data,
      MARGIN = 3,
      STATS = z[, k],
      FUN = "*"
    )
    S_weighted[, , k] <- apply(weighted_data, MARGIN = c(1, 2), FUN = sum) /
      n_K[k]
  }
  S_weighted
  # out <- list(weighted_data=weighted_data, S_k_weighted=S_k_weighted)
  # out
}


log_dens_pro_calculator <- function(data, nu, Sigma, pro, N, K) {

  # Initialize log_dens matrix
  log_dens <- matrix(NA, N, K)

  # Compute log_dens values
  for (k in 1:K) {
    log_dens[, k] <- apply(
      data,
      MARGIN = 3,
      FUN = function(Omega){
        LaplacesDemon::dwishart(
          Omega = Omega,
          nu = nu[k],
          S = Sigma[, , k],
          log = TRUE
        )
      }
    )
  }

  # Compute log_dens_pro values
  log_dens_pro <- sweep(log_dens, 2, log(pro), "+")

  return(log_dens_pro)
}

Q_function_calculator <- function(data, nu, Sigma, pro, z, N, K,LAMBDA) {
  # Calculate log density
  log_dens_pro <- log_dens_pro_calculator(
    data = data,
    nu = nu,
    Sigma = Sigma,
    pro = pro,
    N = N,
    K = K
  )

  penalty_term <- sum(sapply(1:K, function(k)  sum(abs(LAMBDA[,,1]*Sigma[,,k]))))

  # Calculate Q function
  Q_function <- sum(log_dens_pro * z) - penalty_term

  return(Q_function)
}


M_step_sigma_and_nu <- function(data, z, n_K, p, K, nu, LAMBDA, tol, max_iter, N, pro, nu_max) {

  crit_Q_M <- TRUE
  iter_Q_M <- 0
  Q_M <- Q_M_prev <- -.Machine$double.xmax
  # Q_M_vec <- NULL

  # Initialize Sigma and Sigma_inv arrays
  Sigma <- array(0, dim = c(p, p, K))
  Sigma_inv <- array(0, dim = c(p, p, K))


  while (crit_Q_M) {
    # Calculate weighted_S using weighted_sample_S_calculator function
    weighted_S <- weighted_sample_S_calculator(
      data = data,
      z = z,
      n_K = n_K,
      p = p
    )

    # Scale weighted_S by nu
    S_tilde <- sweep(weighted_S, MARGIN = 3, STATS = nu, FUN = "/")

    for (k in 1:K) {

      covglasso_tmp <- covglasso::covglasso(
        S = S_tilde[, , k],
        n = n_K[k], # not used
        lambda = 2 * LAMBDA / (n_K[k] * nu[k]),
        penalize.diag = TRUE
      )

      # Compute Sigma and Sigma_inv for each k
      Sigma[, , k] <- covglasso_tmp$sigma

      Sigma_inv[, , k] <- covglasso_tmp$omega

    # Compute nu for each k
      nu[k] <- tryCatch(
        stats::uniroot(
          f = nu_k_root_finder_equation_cpp,
          interval = c((p - 1) + sqrt(.Machine$double.eps), nu_max),
          data = data,
          p = p,
          z_k = z[, k],
          n_k = n_K[k],
          Sigma_inv_k = Sigma_inv[, , k],
          N = N
        )$root,
        error = function(e)
          print("Error in numerically finding nu_k: values at end points of the function not of opposite sign. Try to increase nu_max")
      )
    }

    # Sanity check: estimated degrees of freedom must be greater or equal than p
    nu <- ifelse(nu < p, p, nu)

    Q_M <- Q_function_calculator(data, nu, Sigma, pro, z,N, K,LAMBDA)

    err_Q_M <- abs(Q_M - Q_M_prev) / (1 + abs(Q_M))
    Q_M_prev <- Q_M
    # Q_M_vec <- c(Q_M_vec, Q_M)

    crit_Q_M <- (err_Q_M > tol & iter_Q_M < max_iter)
    iter_Q_M <- iter_Q_M + 1  # Increment iteration counter
  }

  return(list(Sigma = Sigma, nu = nu))
}


E_step_z_log_lik <- function(data, nu, Sigma, pro, N, K) {

  # Calculate log density with proportionality constant
  log_dens_pro <- log_dens_pro_calculator(
    data = data,
    nu = nu,
    Sigma = Sigma,
    pro = pro,
    N = N,
    K = K
  )

  # Calculate the maximum log density for each row
  zMax <- apply(log_dens_pro, 1, max)

  # Compute the log-likelihood for each obs
  loglik_i <- zMax + log(rowSums(exp(log_dens_pro - zMax)))

  # Update the responsibility matrix
  z <- exp(log_dens_pro - loglik_i)

  return(list(z = z, loglik_i = loglik_i))
}
