compute_bw <- function(Y){
  ## `Y` is a matrix of shape (n, p).
  ## compute bandwidth for each variable (column) of `Y`, as the median of all pairwise difference
  num_sim <- 1000
  bw_list <- lapply(1:ncol(Y), function(idx){
    dist_vec <- rep(NA, num_sim)
    Yj_temp<- stats::na.omit(Y[, idx])
    for (i in 1:num_sim){
      vec <- sample(Yj_temp, 2, replace = FALSE)
      dist_vec[i] <- abs(vec[1] - vec[2])
      }
    bw <- stats::median(dist_vec)  # median of the pairwise distances
    bw <- ifelse(bw > 0, bw, 0.01)
    return (list(bandwidth = bw))
  })

  return (bw_list)
}
