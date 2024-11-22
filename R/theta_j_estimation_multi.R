theta_j_estimation_multi <- function(bj_multi,
                                     B,
                                     cj_multi,
                                     eta_j_multi,
                                     G_multi,
                                     interaction_term,
                                     yj_multi,
                                     adj_col = NULL,
                                     nzero_thres = NULL){
  # multi-sample version theta_j estimation
  # simply concatenate zj and G of all cells by rows, respectively, resulting zj of length Kn and G of dimension (Kn, p^2) or (Kn, p)
  # `bj_multi` is vector of length K (bj for each cell)
  # `cj_multi` is matrix of dimension (n, K) (cj for each cell)
  # `eta_j_multi` is of length K (best eta_j for each cell)
  # `yj_multi` is of dimension (n, K), i.e. all cell series for variable j
  # `G_multi` is a matrix of dimension (Kn, p^2) or (Kn, p). NOTE that `G_multi` is ALREADY CONCATENATED using all cells before passing into this function.

  n <- dim(yj_multi)[1]
  K <- dim(yj_multi)[2]

  zj_multi <- sapply(1:K,
                     function(k){
                       # convert to scalar or vector form
                       Yj <- c(yj_multi[,k])
                       cj <- c(cj_multi[,k])
                       bj <- c(bj_multi[k])
                       eta_j <- c(eta_j_multi[k])

                       zj <- Yj - mean(Yj) - 0.5*n*K*eta_j * cj - B*bj  ### multi-sample version ###
                       return (zj)
                     },
                     simplify = TRUE)  # n x K matrix
  zj <- c(zj_multi)  # concatenate all columns, resulting in a vector of length Kn

  G <- G_multi
  adj_col <- c(adj_col)

  if (any(is.na(G))) {stop("NA found in G")}
  if (any(is.na(zj))) {stop("NA found in zj")}

  if (!is.null(adj_col)){  # if a network is given
    # If a network is given, i.e. which variables affect the current variable j, then implement non-negative linear regression of zj on columns of G corresponding to the interacting variables
    theta_j <- rep(0, ncol(G))
    best_kappa <- 0  # OLS without regularization

    # identify the indices of columns of G that involves the interacting variables
    var_idx <- which(adj_col == 1)  # indices of interacting variables

    if (length(var_idx) == 0) {return (list(best_kappa = best_kappa,
                                            theta_j = theta_j))}  # no interacting variables, return all zeros

    if (interaction_term) {  # G is n x p^2
      idx_mat_tmp <- matrix(1:ncol(G), nrow = sqrt(ncol(G)), ncol = sqrt(ncol(G)), byrow = TRUE)
      col_idx <- idx_mat_tmp[as.matrix(expand.grid(var_idx, var_idx))]  # indices of selected columns of G (all related main and interaction terms)
    } else {  # G is n x p
      col_idx <- var_idx
    }
    G_sub <- G[, col_idx]

    if (length(col_idx) == 1) {G_sub <- cbind(0, G_sub)}  # add pseudo column because glmnet does not support regressing on one single predictor
    if (all(apply(G_sub, 2, stats::var) <= 1e-12)) {
      warning("All columns of G_sub have zero variance, possibly due to cj and theta_j are all zero.",
              "theta_j set to be 0.")
      return (list(best_kappa = best_kappa,
                   theta_j = theta_j))  # cannot perform regression, return all zeros
    }


    # non-negative linear regression of zj on selected G columns
    lm_fit <- glmnet::glmnet(x = G_sub, y = zj, lambda = best_kappa, lower.limits = 0, intercept = FALSE)

    if (length(col_idx) == 1) {
      theta_j[col_idx] <- as.vector(stats::coef(lm_fit)[-c(1,2)])  # exclude the (zero) intercept and pseudo column 1
    } else {
      theta_j[col_idx] <- as.vector(stats::coef(lm_fit)[-1])  # exclude the (zero) intercept
    }
  } else {  # if no network is given, implement non-negative Lasso as described in the paper
    if (all(apply(G, 2, stats::var) <= 1e-12)) {
      warning("All columns of G have zero variance, possibly due to cj and theta_j are all zero.",
              "theta_j set to be 0.")
      return (list(best_kappa = best_kappa,
                   theta_j = rep(0, ncol(G))))  # cannot perform regression, return all zeros
    }

    cv_fit <- glmnet::cv.glmnet(x = G, y = zj, family = "gaussian", alpha = 1, nfolds = min(10, floor(nrow(G)/3)), intercept = FALSE, standardize = TRUE, lower.limits = 0)  # ensure at least 3 obs in a fold

    if (is.null(nzero_thres)) {  # no constraint on the number of nonzero coefficients (i.e. number of selected edges)
      best_kappa <- cv_fit$lambda.1se
      theta_j <- as.vector(stats::coef(cv_fit, s = best_kappa)[-1])  # 1se rule, exclude the (zero) intercept
    } else {  # controlling the number of nonzero coefficients (only used when NO NETWORK IS GIVEN)
      best_kappa <- min(cv_fit$lambda[cv_fit$nzero <= floor(nzero_thres * ncol(G))])
      theta_j <- as.vector(stats::coef(cv_fit, s = best_kappa)[-1])  # exclude the (zero) intercept
    }
  }

  best_kappa <- best_kappa / n  # the regularization parameter is actually n*kappa, refer to Eq (15)
  theta_j[is.na(theta_j)] <- 0  # some coefficients are NA due to duplicated columns of G. At least it happens in lm().

  return (list(best_kappa = best_kappa, theta_j = theta_j))
}

