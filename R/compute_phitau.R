compute_phitau <-
function(dev, hv, devu, hvu, X.rand.disp) {
 g11 <- glm((as.numeric(dev/(1 - hv))) ~ 1, family = Gamma(link=log), weights = as.numeric((1 - hv)/2))
 var.e <- exp(as.numeric(g11$coef[1]))
 g12 <- glm(devu/(1-hvu) ~ X.rand.disp - 1, family = Gamma(link = log), weights = (1-hvu)/2)
 phi <- as.numeric(g12$fitted.value)
 list(var.e = var.e, phi=phi, coef = g12$coef)
}
