.theta_j_estimation <- function(bj,
                                B,
                                cj,
                                eta_j,
                                G,
                                interaction_term,
                                Yj,
                                adj_col = NULL,
                                nzero_thres = NULL){
  n <- length(Yj)

  # convert to scalar or vector form
  Yj <- c(Yj)
  cj <- c(cj)
  bj <- c(bj)
  eta_j <- c(eta_j)

  zj <- Yj - mean(Yj) - 0.5*n*eta_j*cj - B*bj
  adj_col <- c(adj_col)

  if (any(is.na(G))) {stop("NA found in G")}
  if (any(is.na(zj))) {stop("NA found in zj")}

  if (!is.null(adj_col)){  # if a network is given
    # If a network is given, i.e. which variables affect the current variable j, then implement non-negative linear regression of zj on columns of G corresponding to the interacting variables
    theta_j <- rep(0, ncol(G))
    best_kappa <- 0  # OLS without regularization

    # identify the indices of columns of G that correspond to the interacting variables
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

    # non-negative linear regression of zj on selected G columns
    lm_fit <- glmnet::glmnet(x = G_sub, y = zj, lambda = best_kappa, lower.limits = 0, intercept = FALSE)

    if (length(col_idx) == 1) {
      theta_j[col_idx] <- as.vector(stats::coef(lm_fit)[-c(1,2)])  # exclude the (zero) intercept and pseudo column 1
    } else {
      theta_j[col_idx] <- as.vector(stats::coef(lm_fit)[-1])  # exclude the (zero) intercept
    }
  } else {  # if no network is given, implement non-negative Lasso as described in the paper
    cv_fit <- glmnet::cv.glmnet(x = G, y = zj, family = "gaussian", alpha = 1, nfolds = min(10, floor(nrow(G)/3)), intercept = FALSE, standardize = TRUE, lower.limits = 0)  # ensure at least 3 obs in a fold

    if (is.null(nzero_thres)) {  # no constraint on the number of nonzero coefficients (i.e. number of selected edges)
      # select lambda using the 1se rule
      best_kappa <- cv_fit$lambda.1se
      theta_j <- as.vector(stats::coef(cv_fit, s = best_kappa)[-1])  # exclude the (zero) intercept
    } else {  # controlling the number of nonzero coefficients (only used when NO NETWORK IS GIVEN)
      best_kappa <- min(cv_fit$lambda[cv_fit$nzero <= floor(nzero_thres * ncol(G))])
      theta_j <- as.vector(stats::coef(cv_fit, s = best_kappa)[-1])  # exclude the (zero) intercept
    }
  }

  best_kappa <- best_kappa / n  # the regularization parameter is actually n*kappa, refer to Eq (15)
  theta_j[is.na(theta_j)] <- 0  # some coefficients are NA due to duplicated columns of G (at least it happens in `lm()`).

  return (list(best_kappa = best_kappa,
               theta_j = theta_j))
}
