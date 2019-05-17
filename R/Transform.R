Transform <-
function(X, L, var.e, var.u, v, svdVec, svdD, wZt, w) {
    N<-nrow(X)
    A11 <- crossprod(X, X) #Has to be invertible
    A12 <- crossprod(X, L)
    A22 <- crossprod(L) + diag(N)*var.e/var.u
    C22 <- solve(A22-t(A12)%*%solve(A11)%*%A12)*var.e
    rm(list = c("A11", "A22" , "A12"))
    transf <- wZt%*%t( t(svdVec)/sqrt(svdD) ) #Should be Z%*%t( t(svdVec)/sqrt(svdD) ) but corrected below. Saves memory.
    u <- sqrt(w)*transf%*%v #Multiply by sqrt(w) here because wZt is computed as sqrt(w)*t(Z)
    # A vector of ones is given to the following function because wZt is computed as sqrt(w)*t(Z)
    qu <- hat.transf(C22, transf, vc = var.u, k = nrow(wZt), N = nrow(X), w = rep(1, length(w)))
    list(u=u, qu=qu)
}
