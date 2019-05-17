compute_GL <-
function(Z, w) {
    wZt <- sqrt(w)*t(Z)
    G <- crossprod(wZt)
    svdG <- svd(G)
    L <- t(t(svdG$u)*sqrt(svdG$d)) 
    list(L=L, svdVec = svdG$v, svdD = svdG$d, wZt=wZt)
}
