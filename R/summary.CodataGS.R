`summary.CodataGS` <-
  function(object, ...) {
    if (object$Converge == "Not") stop("There is no valid estimate to produce summary statistics")
    #if (is.null(object$fixef)) stop("There in no valid estimate to produce summary statistics")
    #if (is.null(object$SeFe)) stop("There in no valid standard error estimate to produce summary statistics")
    cat("Number of observations in the training set: ", length(object$gEBV), "\n")
    if (!is.null(object$predicted.gEBV)) cat("Number of observations in the test set: ", length(object$predicted.gEBV), "\n")
    cat("Number of markers included: ", length(object$u), "\n")
    cat("\n")
    cat("Estimated fixed effects in the mean part \n")
    cat(round(object$beta, 4), "\n")
    cat("Estimated fixed effects in the model for the random effect variances \n")
    cat("(These estimates are for the linear predictor on a log scale) \n")
    cat(round(object$disp.beta, 4), "\n")
    cat("\n")
    PMSE1 <- mean(object$res.sq/(1-object$hv)^2)
    cat("Mean square error for leave-one-out cross validation \n")
    cat("PMSE(1) = ", PMSE1, "\n")
}

