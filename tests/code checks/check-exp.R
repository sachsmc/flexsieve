
# exp( - (t / 3 + t / 2)) = 1 - p
# log(1 - p) = 2t / 6 + 3t / 6 = 5t/6
# t = -log(1-p) / (5/6)

samdat <- function() {
times <- -log(runif(2000)) / (5/6)
cens <- runif(2000, 0, max(times))

prob1 <- (times / 3) / (times / 3 + times / 2)

ev <- ifelse(rbinom(2000, 1, prob1) == 1, 1, 2)
event <- ifelse(cens < times, 0, ev)

data.frame(ID = c(1:2000, 1:2000), time = c(times, times), status = c(ifelse(event == 1, 1, 0), ifelse(event == 2, 1, 0)),
           trans = c(rep(1, 2000), rep(2, 2000)))

}
oneflex <- function(i){
dataou1 <- list(list(data1 = subset(samdat(), trans == 1), data2 = subset(samdat(), trans == 1)),
                list(data1 = subset(samdat(), trans == 2), data2 = subset(samdat(), trans == 2)))
fittedmods <- NULL
for(i in 1:2) {
    fittedmods[[i]] <- estimate_flexparams(dataou1[[i]], "ID", Surv(time,  status) ~ 1, K = 3, method = "BFGS")
}
nknot <- length(fittedmods[[1]]$knots)
    obj <- function(s){
        cumhaz1 <- Hsurvspline(s, fittedmods[[1]]$optim$par[1:nknot],
                               beta = 0,
                               X = 0,
                               knots = fittedmods[[1]]$knots)

        cumhaz2 <- Hsurvspline(s, fittedmods[[2]]$optim$par[1:nknot],
                               beta = 0,
                               X = 0,
                               knots = fittedmods[[2]]$knots)

        haz1 <- hsurvspline(s, fittedmods[[1]]$optim$par[1:nknot],
                            beta = 0,
                            X = 0,
                            knots = fittedmods[[1]]$knots)

        (exp(-(cumhaz1 + cumhaz2)) * haz1)
    }

    tim <- seq(0.1, 7, by = .1)
    CI.flex <- sapply(tim, function(t) integrate(obj, lower = 0, upper = t)$value)
    CI.flex
}

flexon <- lapply(1:10, oneflex)
    obj2 <- function(s) {
        exp(-(s / 3 + s / 2)) * 1/3
    }
    tim <- seq(0.1, 7, by = .1)
    CI.true <- sapply(tim, function(t) integrate(obj2, 0, t)$value)
    plot(CI.true ~tim, type = "l", ylim = c(0, 1))
    for(j in 1:10) {
    lines(flexon[[j]] ~ tim, col = "blue", lty = 3)
    }
