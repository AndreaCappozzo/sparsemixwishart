# library(abind)
# set.seed(42)
# p <- 10
# Sigma_1 <- diag(p)*runif(n = p,min = 3,max = 7)
# Sigma_2 <- clusterGeneration::genPositiveDefMat(dim = p,covMethod = "eigen")$Sigma
# Sigma_2 <- ifelse(Sigma_2<1e-5,0,Sigma_2)
# Sigma_3 <- clusterGeneration::genPositiveDefMat(dim = p,covMethod = "eigen")$Sigma
#
# n_K <- c(200,400,300)
# class_true <- rep(1:3,n_K)
# z_true <- mclust::unmap(class_true)
# N <- sum(n_K)
# v <- c(20,30,40)
#
# # Sample from Wishart with sparse V param
# X_1 <- rWishart(n = n_K[1], df = v[1], Sigma = Sigma_1)
# X_2 <- rWishart(n = n_K[2], df = v[2], Sigma = Sigma_2)
# X_3 <- rWishart(n = n_K[3], df = v[3], Sigma = Sigma_3)
#
# data <- abind(X_1,X_2,X_3,along = 3)
# K <- 3
# penalty <- 0
# control = sparsemixwishart::EM_controls()
# penalize_diag <- FALSE
# verbose = interactive()
#
