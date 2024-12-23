"Author: Md. Mamunur Rashid
< mamun.stat92@gmail.com > 
MOMLIN R version; June 14, 2024

"

# ===== Define the jet.colors function ========
jet.colors <- function(n) {
  jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                            "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  jet(n)
}


# ====== For L2,1-norm & L1,1-norm ========
updateD <- function(W) {
  if (missing(W)) {
    stop("Error: input W is missing.")
  }
  
  # For L2,1-norm & L1,1-norm
  d <- 1 / sqrt(rowSums(as.matrix(W)^2) + .Machine$double.eps)
  D <- diag(0.5 * d)
  return(D)
}

# Example usage
# W <- matrix(c(1, 2, 3, 4, 5, 6), nrow=3, ncol=2)
# D <- updateD2(W)
# print(D)


# ====== Normalize function ===============
normalize <- function(X, type = 'std') {
  eps <- .Machine$double.eps
  
  if (missing(type)) {
    type <- 'std'
  }
  Y <- switch(type,
              'center' = {
                sweep(X, 2, colMeans(X), '-')
              },
              'median' = {
                sweep(X, 2, matrixStats::colMedians(X), '-')
              },
              'minmax' = {
                X_min <- sweep(X, 2, apply(X, 2, min), '-')
                X_range <- apply(X, 2, max) - apply(X, 2, min) + eps
                sweep(X_min, 2, X_range, '/')
              },
              'norm' = {
                X0 <- sweep(X, 2, colSums(abs(X)) + eps, '/')
                X0_centred <- sweep(X0, 2, colMeans(X0), '-') 
                Y <- sweep(X0_centred, 2, sqrt(colSums(X0_centred^2) + eps), '/') 
                Y
              },
              'std' = {
                X0 <- sweep(X, 2, colSums(abs(X)) + eps, '/')
                scale(X0)
              },
              'none' = {
                X
              },
              stop('Error: Unsupported normalization type.')
  )
  return(Y)
}

# Create a synthetic matrix (N = 10, D = 5)
# X <- matrix(rnorm(50, mean = 10, sd = 5), nrow = 10, ncol = 5)
# rstd_X <- normalize(X, 'rstd')


# ========= Get co-expressed feature by laplacian =========
get_connectivity <- function(Data, alpha = 2) {
  # Data : expression (n * p)
  # alpha : default to 2
  
  # Calculate the correlation matrix or adjancy matrix: A
  A <- cor(Data)^alpha
  A[A <= 0.1]=0  # remove lowest interactions
  # Set the diagonal elements to zero
  diag(A) <- 0
  # Degree matrix
  D <- diag(rowSums(A))
  # Laplacian computation
  L <- D - A 
  return(L)
}

# Example usage
# set.seed(123)  # For reproducibility
# Data <- matrix(rnorm(50, mean = 10, sd = 5), nrow = 10, ncol = 5)
# alpha <- 2
# Penalty <- get_connectivity(Data, alpha)


#  ======== one-hot encode labels ==================
one_hot_encode <- function(labels, num_classes) {
  n <- length(labels)
  L <- matrix(0, nrow = n, ncol = num_classes)
  for (i in 1:n) {
    L[i, labels[i]] <- 1
  }
  L
}


# ======= Class balancing using over-sampling ============
do_oversample <- function(rawData) {
  # rawData is a list with elements X, Y, and Z
  X <- list()
  Y <- list()
  Z <- list()
  
  # Loop through each row of rawData$Y
  for (ii in 1:length(rawData$Y)) {
    X1 <- rawData$X[[ii]]
    Y1 <- rawData$Y[[ii]]
    Z1 <- rawData$Z
    
    tempX <- list()
    tempY <- list()
    tempZ <- list()
    
    for (c in 1:ncol(Z1)) {
      idx_p <- which(Z1[, c] == 1)
      num_p <- length(idx_p)  # Number of positive samples
      idx_n <- which(Z1[, c] == 0)
      num_n <- length(idx_n)  # Number of negative samples
      
      if (num_p > num_n) {
        idx_oversample <- sample(idx_n, num_p - num_n, replace = TRUE)
      } else if (num_n > num_p) {
        idx_oversample <- sample(idx_p, num_n - num_p, replace = TRUE)
      } else {
        idx_oversample <- integer(0)
      }
      
      tempX[[c]] <- rbind(X1, X1[idx_oversample, ])
      tempY[[c]] <- rbind(Y1, Y1[idx_oversample, ])
      tempZ[[c]] <- as.matrix(c(Z1[, c], Z1[idx_oversample, c]))
    }
    
    X[[ii]] <- tempX
    Y[[ii]] <- tempY
    Z[[ii]] <- tempZ
  }
  
  return(list(X = X, Y = Y, Z = Z))
}

# # Example usage
# set.seed(123)  # For reproducibility
# rawData <- trainData
# result <- do_oversample(rawData)


#  ======= Canonical correlations ==============
calcCCC <- function(Data, U, V, scale_type = "none") {
  library(stats)
  
  # Extracting data from the input list
  X <- Data$X
  Y <- Data$Y
  n_class <- Data$n_class
  
  # Pre-allocation of CCCs array (modality x class)
  CCCs <- array(0, dim = c(length(Y), n_class))
  
  for (c in 1:n_class) {
    for (j in 1:length(Y)) {
      
      # Scaling or normalizing data based on the specified scale_type
      if (scale_type == 'none') {
        norm_X <- X[[j]]
        norm_Y <- Y[[j]]
      } else if (scale_type == 'std') {
        norm_X <- normalize(X[[j]], scale_type)
        norm_Y <- normalize(Y[[j]], scale_type)
      } else if (scale_type == 'std1') {
        norm_X <- scale(X[[j]])
        norm_Y <- scale(Y[[j]])
      } else {
        stop("Invalid scale_type specified. Choose 'none', 'std', or 'std1'.")
      }
      # Calculate and store the absolute value of the correlation coefficient
      CCCs[j, c] <- abs(cor(norm_X %*% U[, c], norm_Y %*% V[[j]][, c]))
    }
  }
  return(CCCs)
}


# ============ MOMLIN main function =========
momlin_func <- function(data_train, data_test, opts) {
  library(Matrix)
  
  scale_type = opts$scale_type
  k <- opts$k
  ii <- opts$ii
  
  # weighted Multi-class Sparse Canonical Correlation Analysis
  p <- sapply(data_train$X, ncol)
  q <- sapply(data_train$Y, ncol)
  
  # number of sample
  n = nrow(data_train$X[[1]]) 
  
  # disease/diagnosis class
  n_class <- data_train$n_class
  
  # mRNA weight matrix
  U <- matrix(1, nrow = p[1], ncol = n_class)
  
  # miRNA, Methy, & CNV weights matrix
  V <- lapply(q, function(dim) matrix(1, nrow = dim, ncol = n_class))
  
  # Z weights
  W <- matrix(1, nrow = 1, ncol = n_class)
  
  # Calculate connectivity penalty
  Lgn_u <- lapply(data_train$X, function(X) get_connectivity(X, 2))
  
  # class balancing by oversampling
  oversampled <- do_oversample(data_train)
  X <- oversampled$X
  Y <- oversampled$Y
  Z <- oversampled$Z
  
  # initial covariance & scale matrix
  # Initialize nested lists
  XX <- lapply(1:length(Y), function(i) vector("list", n_class))
  XY <- lapply(1:length(Y), function(i) vector("list", n_class))
  YY <- lapply(1:length(Y), function(i) vector("list", n_class))
  YX <- lapply(1:length(Y), function(i) vector("list", n_class))
  ZZ <- lapply(1:length(Y), function(i) vector("list", n_class))
  XZ <- lapply(1:length(Y), function(i) vector("list", n_class))
  YZ <- lapply(1:length(Y), function(i) vector("list", n_class))
  ZX <- lapply(1:length(Y), function(i) vector("list", n_class))
  ZY <- lapply(1:length(Y), function(i) vector("list", n_class))
  
  # weighting function
  sigma_xyk <- lapply(1:length(Y), function(i) vector("list", n_class))
  sigma_xz <- lapply(1:length(Y), function(i) vector("list", n_class))
  sigma_yzk <- lapply(1:length(Y), function(i) vector("list", n_class))
  
  
  for (c in seq_len(n_class)) {
    for (i in seq_along(Y)) {
      # data normalization
      X[[i]][[c]] <- normalize(X[[i]][[c]], scale_type)
      Y[[i]][[c]] <- normalize(Y[[i]][[c]], scale_type)
      
      # pre-compuations
      XX[[i]][[c]] <- crossprod(X[[i]][[c]])
      XY[[i]][[c]] <- crossprod(X[[i]][[c]], Y[[i]][[c]])
      YY[[i]][[c]] <- crossprod(Y[[i]][[c]])
      YX[[i]][[c]] <- t(XY[[i]][[c]])
      
      ZZ[[i]][[c]] <- crossprod(Z[[i]][[c]])
      XZ[[i]][[c]] <- crossprod(X[[i]][[c]], Z[[i]][[c]])
      YZ[[i]][[c]] <- crossprod(Y[[i]][[c]], Z[[i]][[c]])
      ZX[[i]][[c]] <- t(XZ[[i]][[c]])
      ZY[[i]][[c]] <- t(YZ[[i]][[c]])
      
      # scale weights
      U[, c] <- U[, c] / norm(X[[i]][[c]] %*% U[, c],"2") 
      V[[i]][, c] <- V[[i]][, c] / norm(Y[[i]][[c]] %*% V[[i]][, c],"2")
      W[, c] <- W[, c] / norm(Z[[i]][[c]] %*% W[, c],"2")
      
      Xu <- X[[i]][[c]] %*% U[, c]
      Yv <- Y[[i]][[c]] %*% V[[i]][, c]
      Zw <- Z[[i]][[c]] %*% W[, c]
      
      sigma_xyk[[i]][[c]] <- 1 / norm(Xu - Yv, "2")
      sigma_xz[[i]][[c]] <- 1 / norm(Xu - Zw, "2")
      sigma_yzk[[i]][[c]] <- 1 / norm(Yv - Zw, "2")
    }
  }
  
  
  # Parameters
  alpha_u <- opts$alpha_u
  alpha_v <- list(opts$alpha_v1, opts$alpha_v2, opts$alpha_v3, opts$alpha_v4)
  beta <- opts$beta
  
  max_Iter <- 50
  iter <- 0
  tol <- 1e-5
  tu <- Inf
  tv <- Inf
  
  err <- 0.005
  
  objFun_train <- objFun_test <- numeric(max_Iter)
  
  diff_obj_train <- numeric()
  diff_obj_test <- numeric()
  loss <- list()
  
  while (iter < max_Iter && (tu > tol || tv > tol)) {
    iter <- iter + 1
    
    U_old <- U
    V_old <- V
    
    
    # ==== Update U =====
    for (c in seq_len(n_class)) {
      Du_c <- updateD(U[, c])
      
      for (i in seq_along(Y)) {
        Fu <- alpha_u * beta * Du_c + alpha_u * (1 - beta) * Lgn_u[[i]] + XX[[i]][[c]]
        bu <- sigma_xyk[[i]][[c]] * XY[[i]][[c]] %*% V[[i]][, c] + sigma_xz[[i]][[c]] * XZ[[i]][[c]]
        
        U[, c] <- solve(Fu, bu)  # (Fu)^-1*bu 
        # Scale each u_c
        U[, c] <- U[, c] / norm(X[[i]][[c]] %*% U[, c],'2')
      }
    }
    
    
    # ===== Update V_k ========
    for (c in seq_len(n_class)) {
      
      for (i in seq_along(Y)) {
        Dv_c <- updateD(V[[i]][, c])
        
        Fv <- alpha_v[[i]] * Dv_c + YY[[i]][[c]]
        bv <- sigma_xyk[[i]][[c]] * YX[[i]][[c]] %*% U[, c] + sigma_yzk[[i]][[c]] * YZ[[i]][[c]]
        
        V[[i]][, c] <- solve(Fv, bv)   # (Fv)^-1*bv 
        # Scale each V.k_c
        V[[i]][, c] <- V[[i]][, c] / norm(Y[[i]][[c]] %*% V[[i]][, c], '2')
      }
    }
    
    
    # ===== Update W ==========
    for (c in seq_len(n_class)) {
      for (i in seq_along(Y)) {
        Fw <- ZZ[[i]][[c]]
        # bw <- ZX[[i]][[c]] %*% U[, c] + ZY[[i]][[c]] %*% V[[i]][, c]
        bw <- sigma_xz[[i]][[c]] * ZX[[i]][[c]] %*% U[, c] + sigma_yzk[[i]][[c]] * ZY[[i]][[c]] %*% V[[i]][, c]
        
        W[, c] <- solve(Fw, bw)
        W[, c] <- W[, c] / norm(Z[[i]][[c]] %*% W[, c], '2')
      }
    }
    
    
    # Update the dynamics weights
    for (c in seq_len(n_class)) {
      for (i in seq_along(Y)) {
        Xu <- X[[i]][[c]] %*% U[, c]
        Yv <- Y[[i]][[c]] %*% V[[i]][, c]
        Zw <- Z[[i]][[c]] %*% W[, c]
        
        sigma_xyk[[i]][[c]] <- 1 / norm(Xu - Yv, '2')
        sigma_xz[[i]][[c]] <- 1 / norm(Xu - Zw, '2')
        sigma_yzk[[i]][[c]] <- 1 / norm(Yv - Zw, '2')
      }
    }
    
    
    # Iteration termination
    tu <- max(abs(U - U_old))
    tv <- sapply(seq_along(Y), function(i) max(abs(V[[i]] - V_old[[i]])))
    tv <- max(tv)
    # ========
    tuv = max(tu, tv)
    
    # Cost function and check convergence
    
    for (i in seq_along(Y)) {
      objFun_train[iter] <- -sigma_xyk[[i]][[1]] * sum(diag(crossprod(U, crossprod(data_train$X[[i]], data_train$Y[[i]] %*% V[[i]])))) -
                                                      sigma_xz[[i]][[1]] * sum(diag(crossprod(U, crossprod(data_train$X[[i]], data_train$Z)))) -
                                                      sigma_yzk[[i]][[1]] * sum(diag(crossprod(V[[i]], crossprod(data_train$Y[[i]], data_train$Z)))) +
                                                      alpha_u * beta * sum(sqrt(rowSums(U^2))) +
                                                      alpha_u * (1 - beta) * sum(diag(crossprod(U, Lgn_u[[i]] %*% U))) +
                                                      alpha_v[[i]] * sum(sqrt(rowSums(V[[i]]^2)))
      
      objFun_test[iter] <- -sigma_xyk[[i]][[1]] * sum(diag(crossprod(U, crossprod(data_test$X[[i]], data_test$Y[[i]] %*% V[[i]])))) -
                                                     sigma_xz[[i]][[1]] * sum(diag(crossprod(U, crossprod(data_test$X[[i]], data_test$Z)))) -
                                                     sigma_yzk[[i]][[1]] * sum(diag(crossprod(V[[i]], crossprod(data_test$Y[[i]], data_test$Z)))) +
                                                     alpha_u * beta * sum(sqrt(rowSums(U^2))) +
                                                     alpha_u * (1 - beta) * sum(diag(crossprod(U, Lgn_u[[i]] %*% U))) +
                                                     alpha_v[[i]] * sum(sqrt(rowSums(V[[i]]^2)))
    }
    
    
    if (iter != 1) {
      if (ii == 1) {
        if (k == 1) {
          diff_obj_train[iter] <- abs((objFun_train[iter] - objFun_train[iter - 1]) / n)
          diff_obj_test[iter] <- abs((objFun_test[iter] - objFun_test[iter - 1]) / n)
          
          loss$train <- diff_obj_train
          loss$test <- diff_obj_test
          loss$iter <- seq_len(iter)
        }
      }
    }
  }
  
  list(U = U, V = V, loss = loss)
}

# End
