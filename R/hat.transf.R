hat.transf <-
function(C22, transf, vc, k, N, w) {
	qu <- numeric(k)
	middle <- diag(N) - C22/vc
	svdmid <- svd(middle)
	M <- t(svdmid$u %*% diag(sqrt(svdmid$d)))
	qu <- 1 - colSums(tcrossprod(M, transf)**2)/w
	#qu <- 1 - colSums(tcrossprod(M, transf)**2)
	return(qu)
}
