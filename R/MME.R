MME <-
function(y, X, Z, var.e, var.u) {
 n <- length(y)
 k <- ncol(Z)
 AugXZ <- cbind(X,Z) #Changed from cBind in version 1.1
 XX1 <- Matrix(0, nrow = ncol(Z), ncol = ncol(X))
 ZZ2 <- Diagonal(ncol(Z))
 ZZ2 <- cbind(XX1, ZZ2) #Changed from cBind in version 1.1
 AugXZ <- rbind(AugXZ, ZZ2) #Changed from rBind in version 1.1
 rm(list = c("XX1", "ZZ2"))
 Augy <- c(y, rep(0, ncol(Z)))
 w <- c(rep(1/sqrt(var.e), n), rep(1/sqrt(var.u), k))		
 SQR <- qr(AugXZ*w)
 est <- as.numeric(qr.coef(SQR, y = Augy*w)) 
 hv <- rowSums(qr.qy(SQR, diag(1, nrow = length(Augy), ncol = ncol(AugXZ)))^2)
 resid <- y - X%*%est[1:ncol(X)] - Z%*%est[-c(1:ncol(X))]
 list(beta=est[1:ncol(X)], v=est[-c(1:ncol(X))], hv=hv[1:n], dev = resid^2)
}
