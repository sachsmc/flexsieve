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


b1 <- matrix(c(1, 0, 0, .5, rep(0, 4)), nrow = 8)
b2 <- matrix(c(1, 0, 0, 0, rep(0, 4)), nrow = 8)
g1 <- 3
g2 <- 3
n <- 5e2
W<-rnorm(n,2,0.5)
S1<-rnorm(n,0.3+W,0.75)
S1 <- S1 - mean(S1)
tr<-rbinom(n,1,0.5)
X3 <- matrix(rnorm(n * 3, mean = 0, sd = .1), ncol = 3)
W<-rnorm(n,2,0.5)

matx<-cbind(1,tr,S1,tr*S1,W,X3)
Touts<-sampletimesfunc(matx,b1,g1,b2,g2)
hist(Touts)


s.obs <- sort(S1)

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
plot(drd.s ~ s.obs, type = "l")
plot((ci.1.trt - ci.1.pla) ~ s.obs, type = "l")
lines((ci.2.trt - ci.2.pla) ~ s.obs, type = "l")

drd <- mean(drd.s)

rd.curve <- as.vector(1 - St(2, Z.trt, b1, g1, b2, g2) ) -
    as.vector(1 - St(2, Z.pla, b1, g1, b2, g2) )

plot(rd.curve ~ s.obs, type = "l")
