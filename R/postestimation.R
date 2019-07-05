#' Compute the cumulative incidence
#'
#' @param estbetas Fitted flexsurv object as returned by \link{estimate_flexparams}
#' @param timet The time at which the cumulative incidence is to be computed
#' @param datain data frame containing all variables in the model
#'
#' @export

cumincfun <- function(estbetas, timet, datain) {

    pred.form <- update(estbetas$formula, NULL ~ .)
    datamat <- model.matrix(pred.form, datain)

    apply(datamat, MAR = 1, FUN = function(x) {

        obj <- function(s){
            haz <- hsurvspline(s, estbetas$optim$par[1:length(estbetas$knots)],
                               beta = estbetas$optim$par[-c(1:length(estbetas$knots))],
                               X = x,
                               knots = estbetas$knots)
            cumhaz <- Hsurvspline(s, estbetas$optim$par[1:length(estbetas$knots)],
                                  beta = estbetas$optim$par[-c(1:length(estbetas$knots))],
                                  X = x,
                                  knots = estbetas$knots)
            haz * exp(-cumhaz)
        }

        integrate(obj, lower = 0, upper = timet)$value

    })

}

#' G compute the cumulative incidence function
#'
#' assuming S(1), W are jointly normal
#' sample from W | S(1), estimate CI, then average over the Ws
#' at values of S(1) that are observed in the vaccine arm
#'
#' @param dataou1 Augmented data as returned by \link{augment_data}
#' @param fittedmod Parameter estimates as returned by \link{estimate_flexparams}
#' @param timet numeric time at which cumulative incidence is computed
#' @param treatment.var string denoting the name of the treatment variable
#' @param surrogate.var string denoting the name of the surrogate variable
#' @param nn number of values to sample for each observed S.1
#'
#' @return A data frame with the variables that go into the model, the estimated cumulative
#' incidence and an indicator for the group to average over
#'
#'
#' @export


gcompute_cuminc <- function(dataou1, fittedmods, timet, treatment.var, surrogate.var, nn = 100) {

    w.vars <- setdiff(all.vars(update(fittedmods[[1]]$formula, NULL ~ .)), c(treatment.var, surrogate.var))

    normdat <- rbind(dataou1[[1]]$data1[, c(surrogate.var, w.vars)], dataou1[[2]]$data1[, c(surrogate.var, w.vars)])
    mu.norm <- sapply(normdat, mean)
    Sigma.norm <- cov(normdat)

    at.S <- normdat[, surrogate.var]

    ## sample W | S(1) for each S(1) in at.S
    w1.dex <- c(1:ncol(normdat))[-1]
    mu.cond <- t(sapply(at.S, function(ss) mu.norm[w1.dex] + Sigma.norm[w1.dex, 1] * (ss - mu.norm[1]) / (Sigma.norm[1, 1])))
    Sig.cond <- Sigma.norm[w1.dex, w1.dex] - Sigma.norm[w1.dex, 1] %*% solve(Sigma.norm[1, 1]) %*% Sigma.norm[1, w1.dex]

    W.samp <- do.call(rbind, lapply(1:nrow(mu.cond), function(x){
        dubus <- as.data.frame(rmvnorm(nn, mean = mu.cond[x, ], sigma = Sig.cond))
        dubus <- rbind(dubus, dubus)
        dubus[[surrogate.var]] <- at.S[x]
        dubus[[treatment.var]] <- rep(c(0, 1), each = nn)
        dubus$group <- x
        dubus
    }))

    stand.X <- model.matrix(update(fittedmods[[1]]$formula, NULL ~ .), W.samp)
    nknot <- length(fittedmods[[1]]$knots)
    nknot2 <- length(fittedmods[[2]]$knots)
    CI.cause1 <- function(i) {
        obj <- function(s){
            cumhaz1 <- Hsurvspline(s, fittedmods[[1]]$optim$par[1:nknot],
                                   beta = fittedmods[[1]]$optim$par[-c(1:nknot)],
                                   X = stand.X[i, , drop = FALSE],
                                   knots = fittedmods[[1]]$knots)

            cumhaz2 <- Hsurvspline(s, fittedmods[[2]]$optim$par[1:nknot2],
                                   beta = fittedmods[[2]]$optim$par[-c(1:nknot2)],
                                   X = stand.X[i, , drop = FALSE],
                                   knots = fittedmods[[2]]$knots)

            haz1 <- hsurvspline(s, fittedmods[[1]]$optim$par[1:nknot],
                                beta = fittedmods[[1]]$optim$par[-c(1:nknot)],
                                X = stand.X[i, , drop = FALSE],
                                knots = fittedmods[[1]]$knots)

            (exp(-(cumhaz1 + cumhaz2)) * haz1)
        }

        obj2 <- function(s){
            cumhaz1 <- Hsurvspline(s, fittedmods[[1]]$optim$par[1:nknot],
                                   beta = fittedmods[[1]]$optim$par[-c(1:nknot)],
                                   X = stand.X[i, , drop = FALSE],
                                   knots = fittedmods[[1]]$knots)

            cumhaz2 <- Hsurvspline(s, fittedmods[[2]]$optim$par[1:nknot2],
                                   beta = fittedmods[[2]]$optim$par[-c(1:nknot2)],
                                   X = stand.X[i, , drop = FALSE],
                                   knots = fittedmods[[2]]$knots)

            haz2 <- hsurvspline(s, fittedmods[[2]]$optim$par[1:nknot2],
                                beta = fittedmods[[2]]$optim$par[-c(1:nknot2)],
                                X = stand.X[i, , drop = FALSE],
                                knots = fittedmods[[2]]$knots)

            (exp(-(cumhaz1 + cumhaz2)) * haz2)
        }

        obj3 <- function(s){
            cumhaz1 <- Hsurvspline(s, fittedmods[[1]]$optim$par[1:nknot],
                                   beta = fittedmods[[1]]$optim$par[-c(1:nknot)],
                                   X = stand.X[i, , drop = FALSE],
                                   knots = fittedmods[[1]]$knots)

            cumhaz2 <- Hsurvspline(s, fittedmods[[2]]$optim$par[1:nknot2],
                                   beta = fittedmods[[2]]$optim$par[-c(1:nknot2)],
                                   X = stand.X[i, , drop = FALSE],
                                   knots = fittedmods[[2]]$knots)

            (exp(-(cumhaz1 + cumhaz2)))
        }

        c(integrate(obj, lower = 0, upper = timet)$value,
          integrate(obj2, lower = 0, upper = timet)$value,
          obj3(timet))

    }

    cuminces <- sapply(1:nrow(stand.X), CI.cause1)
    W.samp$cuminc.cause1 <- cuminces[1, ]
    W.samp$cuminc.cause2 <- cuminces[2, ]
    W.samp$surv.overall <- cuminces[3, ]
    W.samp


}
