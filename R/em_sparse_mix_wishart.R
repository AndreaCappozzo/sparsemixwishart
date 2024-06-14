em_sparse_mix_wishart <- function(data,
                                  K,
                                  penalty,
                                  penalize_diag,
                                  control,
                                  data_dim,
                                  hc_init) {

  p <- data_dim[1]
  N <- data_dim[3]

  # if (!is.matrix(penalty)) penalty <- matrix(penalty, nrow = p, ncol = p) # FIXME need to investigate how to give a matrix to covglasso in order to have different entry-wise penalties
  # if (penalize_diag == FALSE) diag(penalty) <- 0

  # store EM parameters
  tol <- control$tol
  max_iter <- control$max_iter
  n_random_start <-  control$n_random_start

  # initialization of z -----------------------------------------------------------------
  class_init <- if(is.null(hc_init)) {
    sample(1:K, size = N, replace = TRUE)
  } else {
    stats::cutree(hc_init, k = K)
  }
  z <- mclust::unmap(class_init)
  nu <- rep(p,K)
  n_K <- colSums(z)
  weighted_S <- weighted_sample_S_calculator(data = data,z = z,n_K = n_K,p = p)
  S_tilde <- sweep(weighted_S, MARGIN = 3,STATS = nu,FUN = "/")
  Sigma <- Sigma_inv <- array(dim = c(p,p,K))

  for (k in 1:K) {
    Sigma[, , k] <- covglasso::covglasso(S = S_tilde[, , k],
                                           n = n_K[k], # not used
                                           lambda = 2 * penalty / (n_K[k] * nu[k]))$sigma
  }


  iter <- 0
  loglik_pen <- loglik_pen_prev <- -.Machine$integer.max / 2
  loglik_pen_vec <- NULL
  crit <- TRUE

  # ME algorithm

  while (crit) {

    # M step ------------------------------------------------------------------

    #tau
    n_K <- colSums(z)
    pro <- n_K / N

    # Sigma and nu
    m_step_output_Sigma_nu <- M_step_sigma_and_nu(data, z, n_K, p, K, nu, penalty, tol, max_iter, N, pro)

    Sigma <- m_step_output_Sigma_nu$Sigma
    nu <- m_step_output_Sigma_nu$nu

    # E step ------------------------------------------------------------------
    e_step_output <- E_step_z_log_lik(data, nu, Sigma, pro, N, K)
    z <- e_step_output$z

    # Log lik and objective function calculation ------------------------------
    loglik <- sum(e_step_output$loglik_i)

    penalty_term <- penalty * sum(sapply(1:K, function(k)  sum(abs(Sigma[,,k]))))

    loglik_pen <- loglik - penalty_term

    # Check convergence -------------------------------------------------------

    err <- abs(loglik_pen - loglik_pen_prev) / (1 + abs(loglik_pen))
    loglik_pen_prev <- loglik_pen
    loglik_pen_vec <- c(loglik_pen_vec, loglik_pen)
    iter <- iter + 1
    crit <- (err > tol & iter < max_iter)

  }

  # Results collection

  parameters <- list(pro=pro, Sigma=Sigma, nu=nu)
  n_par_pro <- K - 1
  n_par_nu <- K
  n_par_Sigma <- p * K + sum(apply(Sigma,
                                   3, function(A) A[upper.tri(A)] != 0))
  bic_final <- 2 * loglik - (n_par_pro + n_par_nu + n_par_Sigma) * log(N)

  OUT <-
    list(
      loglik = loglik,
      loglik_pen = loglik_pen,
      parameters = parameters,
      z = z,
      K=K,
      classification = mclust::map(z),
      bic = bic_final,
      n_params = c(pro = n_par_pro, nu = n_par_nu, Sigma = n_par_Sigma),
      penalty_term = penalty_term,
      obj_func_trace = loglik_pen_vec,
      iter = iter
    )

  return( OUT )

}
