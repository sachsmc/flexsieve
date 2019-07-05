Cuminc1 <- function(t, Z, beta1, gamma1, beta2, gamma2) {

    obj <- function(u) {
        sapply(u, function(uu) St(uu, Z, beta1, gamma1, beta2, gamma2) * haz11(uu, Z, beta1, gamma1))
    }
    integrate(obj, 0, t)$value

}

Cuminc2 <- function(t, Z, beta1, gamma1, beta2, gamma2) {

    obj <- function(u) {
        sapply(u, function(uu) St(uu, Z, beta1, gamma1, beta2, gamma2) * haz11(uu, Z, beta2, gamma2))
    }
    integrate(obj, 0, t)$value

}


true.curves <- function(b1, b2, g1, g2, n = 5e2, scenario) {

W<-rnorm(n,2,0.5)
S1<-rnorm(n,0.3+W,0.75)
S1 <- S1 - mean(S1)
tr<-rbinom(n,1,0.5)
X3 <- matrix(rnorm(n * 3, mean = 0, sd = .1), ncol = 3)
W<-rnorm(n,2,0.5)

matx<-cbind(1,tr,S1,tr*S1,X3, W)
Touts<-sampletimesfunc(matx,b1,g1,b2,g2)

s.obs <- seq(-3, 3, by = .01)
############## intercept, trt, s, s * trt, X1, X2, X3, W
Z.trt <- cbind(1, 1, s.obs, s.obs, 0, 0, 0, 0)
Z.pla <- cbind(1, 0, s.obs, 0, 0, 0, 0, 0)

ci.1.trt <- sapply(1:nrow(Z.trt), function(i) {
    Cuminc1(2, Z.trt[i, ], b1, g1, b2, g2)
})
ci.1.pla <- sapply(1:nrow(Z.pla), function(i) {
    Cuminc1(2, Z.pla[i, ], b1, g1, b2, g2)
})
ci.2.trt <- sapply(1:nrow(Z.trt), function(i) {
    Cuminc2(2, Z.trt[i, ], b1, g1, b2, g2)
})
ci.2.pla <- sapply(1:nrow(Z.pla), function(i) {
    Cuminc2(2, Z.pla[i, ], b1, g1, b2, g2)
})

drd.s <- (ci.1.trt - ci.1.pla) - (ci.2.trt - ci.2.pla)

drd <- mean(drd.s)

rd.curve <- as.vector(1 - St(2, Z.trt, b1, g1, b2, g2) ) -
    as.vector(1 - St(2, Z.pla, b1, g1, b2, g2) )



list(s = s.obs, drd.s = drd.s, rd.curve = rd.curve, drd = drd)

}
