genomicEBV.w.codata <-
function(y, X, Z, X.SNPcodata, Z.test = NULL, max.iter=100, conv.crit=1e-5) {
  tol.err=1e-6
  #Scale Z according to van Raden 2008
  freq1 <- colMeans(Z)/2
  Z <- scaleZ(Z, freq1)
  if ( !is.null(Z.test) ) Z.test <- scaleZ(Z.test, freq1)
  #Set starting values. Same as hglm package
  g1 <- lm(y ~ X - 1) 
  var.u0 <- (var.e <- as.numeric(.6*summary(g1)$sigma))*.66
  w <- rep(var.u0, ncol(Z))
  w1 <- rep(1, ncol(Z))
  u0 <- rep(0, ncol(Z))
  conv = 999
  iter = 0
  while ( conv > conv.crit & iter < max.iter) {
    iter = iter + 1
    GL <- compute_GL(Z=Z, w=w) 
    mme.obj <- MME(y = y, X = X, Z=GL$L, var.e = var.e, var.u = 1) 
    u.transf <- Transform(X = X, L = GL$L, var.e = var.e, var.u =  1, v = mme.obj$v, svdVec = GL$svdVec, svdD = GL$svdD, wZt = GL$wZt, w = w)
    devu = u.transf$u^2 #Simplified from version 1.1
    phitau <- compute_phitau(dev = mme.obj$dev , hv = mme.obj$hv , devu = devu , hvu = u.transf$qu, X.rand.disp = X.SNPcodata)
    ###
    var.e <- phitau$var.e
    w <- phitau$phi
    w[w < tol.err] <- tol.err
    ###
    conv <- var( abs( u.transf$u - u0 ))
    cat("Iteration", iter, ". Convergence value", as.numeric(conv),"\n")
    u0 <- u.transf$u
  }
  #############
  predicted.gEBV <- NULL
  if(conv < conv.crit) {
    Converge="converged"
    gEBV <- Z %*% u0 #Scaled matrix Z used here
    if ( !is.null(Z.test) ) predicted.gEBV <- Z.test %*% u0
    }
  else {
    Converge="Not"
    gEBV <- NULL
  }
  #Simplified output from version 1.1 with u = u0
  #Class added from version 1.3
  # hv and res.sq added to the list from version 1.3
  val <- list(gEBV = gEBV, predicted.gEBV = predicted.gEBV, w = w, u = u0, beta = mme.obj$beta, disp.beta = phitau$coef, Converge = Converge, iter = iter, hv = mme.obj$hv, res.sq = mme.obj$dev)
  class(val) <- "CodataGS"
  return(val)
}
