dex.ica.R.par = function(X, n.comp, tol, fun, alpha, maxit, verbose, w.init) {
  min.tol = 1
  Diag <- function(d) if (length(d) > 1L)
    diag(d)
  else as.matrix(d)
  p <- ncol(X)
  W <- w.init
  sW <- La.svd(W)
  W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W
  W1 <- W
  lim <- rep(1000, maxit)
  it <- 1
  if (fun == "logcosh") {
    if (verbose)
      message("Symmetric FastICA using logcosh approx. to neg-entropy function")
    while (lim[it] > tol && it < maxit) {
      wx <- W %*% X
      gwx <- tanh(alpha * wx)
      v1 <- gwx %*% t(X)/p
      g.wx <- alpha * (1 - (gwx)^2)
      v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
      W1 <- v1 - v2
      sW1 <- La.svd(W1)
      W1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% W1
      lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
      if(lim[it + 1] < min.tol) min.tol = lim[it + 1]
      W <- W1
      if (verbose)
        message("Iteration ", it, " tol = ", format(lim[it + 1]))
      it <- it + 1
    }
  }
  if (fun == "exp") {
    if (verbose)
      message("Symmetric FastICA using exponential approx. to neg-entropy function")
    while (lim[it] > tol && it < maxit) {
      wx <- W %*% X
      gwx <- wx * exp(-(wx^2)/2)
      v1 <- gwx %*% t(X)/p
      g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
      v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
      W1 <- v1 - v2
      sW1 <- La.svd(W1)
      W1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% W1
      lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
      if(lim[it + 1] < min.tol) min.tol = lim[it + 1]
      W <- W1
      if (verbose)
        message("Iteration ", it, " tol = ", format(lim[it + 1]))
      it <- it + 1
    }
  }
  return(list("W" = W, "iterations" = it, "actual.tol" = lim[it], "min.obs.tol" = min.tol))
}
