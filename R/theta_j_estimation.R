theta_j_estimation <- function(bj,
                               B,
                               cj,
                               eta_j,
                               G,
                               interaction_term,
                               Yj,
                               adj_col = NULL,
                               nzero_thres = NULL){
  ## TODO: consider changing `adj_matrix` to p^2 x p for interaction_terms = TRUE so that no hierarchical structure is assumed
  ## TODO: OR consider using lasso on the selected columns of G when incorporating adjacency matrix into algorithm with interaction terms
  ## TODO: consider only regress on half of the columns of G to reduce multi-colinearity and improve the results


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

      ## model selection alternative 1: additional selection by lm()
      ## use the lambda such that the positive theta_j elements from Lasso are also all significant in an additional lm() call
      ## way 1 (focus on theta_j)
      # G_sub <- G[,theta_j > 0]
      # lm_fit <- lm(zj ~ -1 + G_sub)
      # selector <- (summary(lm_fit)$coefficients[,"Pr(>|t|)"] <= 0.01)
      # theta_j[theta_j > 0][!selector] <- 0
      # ## way 2 (focus on lambda)
      # lambda_candidates <- rep(NA, length(cv_fit$lambda))
      # for (q in 1:length(cv_fit$lambda)){
      #   lambda_temp <- cv_fit$lambda[q]
      #   theta_j_temp <- as.vector(stats::coef(cv_fit, s = lambda_temp)[-1])  # exclude the (zero) intercept
      #   if (sum(theta_j_temp > 0) == 0) {next}  # if no element in theta_j is selected, skip to next lambda value
      #   G_sub_temp <- G[, theta_j_temp > 0]
      #   lm_fit_temp <- lm(zj ~ -1 + G_sub_temp)
      #   if (all(summary(lm_fit_temp)$coefficients[,"Pr(>|t|)"] <= 0.05)){
      #     # if all positive elements of theta_j are still significant in the lm() of zj on columns of G corresponding to them, then this lambda_temp is valid
      #     lambda_candidates[q] <- lambda_temp
      #   }
      # }
      # lambda_candidates <- na.omit(lambda_candidates)
      # if (length(lambda_candidates) == 0){lambda_candidates <- cv_fit$lambda.1se}  # if no lambda satisfies this criterion, use 1se rule
      # lambda_candidates <- lambda_candidates[lambda_candidates >= cv_fit$lambda.1se]  # no smaller than the one selected by 1se rule
      # best_kappa <- min(lambda_candidates)
      # theta_j <- as.vector(stats::coef(cv_fit, s = best_kappa)[-1])  # exclude the (zero) intercept

      ## model selection alternative 2: compute BIC  (https://doi.org/10.1214/009053607000000127)
      # lm_fit <- lm(zj ~ -1 + G)
      # sigma_sq <- (summary(lm_fit)$sigma)^2
      # BICs <- cv_fit$cvm/sigma_sq + log(n)/n * cv_fit$nzero
      # best_kappa <- cv_fit$lambda[which.min(BICs)]
      # theta_j <- as.vector(stats::coef(cv_fit, s = best_kappa)[-1])  # exclude the (zero) intercept

    } else {  # controlling the number of nonzero coefficients (only used when NO NETWORK IS GIVEN)
      best_kappa <- min(cv_fit$lambda[cv_fit$nzero <= floor(nzero_thres * ncol(G))])
      theta_j <- as.vector(stats::coef(cv_fit, s = best_kappa)[-1])  # exclude the (zero) intercept
    }
  }

  best_kappa <- best_kappa / n  # the regularization parameter is actually n*kappa, refer to Eq (15)
  theta_j[is.na(theta_j)] <- 0  # some coefficients are NA due to duplicated columns of G. At least it happens in lm().

  return (list(best_kappa = best_kappa,
               theta_j = theta_j))
}
