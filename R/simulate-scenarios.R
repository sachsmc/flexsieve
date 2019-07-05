switch_scenario <- function(scenario) {


    ## NULL scenario
    if(scenario == "0") {
        b1 <- matrix(c(1, 1, 0, 0, rep(0, 4)), nrow = 8)
        b2 <- matrix(c(1, 1, 0, 0, rep(0, 4)), nrow = 8)
        g1 <- 3
        g2 <- 3
    } else if (scenario == "1") {
        ## Scenario 1: drd 0 on average, but varies with S.1

        b1 <- matrix(c(1.3, 0, 0, .5, rep(0, 4)), nrow = 8)
        b2 <- matrix(c(1.3, 0, 0, -.5, rep(0, 4)), nrow = 8)
        g1 <- 3
        g2 <- 3
    } else if (scenario == "2") {
        ## Scenario 2:
        b1 <- matrix(c(1, 1, 0, 0, rep(0, 4)), nrow = 8)
        b2 <- matrix(c(1, 0, 0, 0, rep(0, 4)), nrow = 8)
        g1 <- 3
        g2 <- 3
    } else if (scenario == "3") {

        ## Scenario 3:

        b1 <- matrix(c(1, 0, 0, .5, rep(0, 4)), nrow = 8)
        b2 <- matrix(c(1, 0, 0, 0, rep(0, 4)), nrow = 8)
        g1 <- 3
        g2 <- 3
    } else if (scenario == "4") {

        b1 <- matrix(c(1, 1, 0, 0, rep(0, 4)), nrow = 8)
        b2 <- matrix(c(1, 1, 0, 0, rep(0, 4)), nrow = 8)
        g1 <- 3
        g2 <- 3
    } else if (scenario == "5") {

        b1 <- matrix(c(1, 0, 0, 0.5, 1, 1, 0, 0), nrow = 8)
        b2 <- matrix(c(1, 0, 0, 0, 1, 1, 0, 0), nrow = 8)
        g1 <- 3
        g2 <- 3

    } else if (scenario == "6") {

        b1 <- matrix(c(1, 0, 0, 0.5, 1, 1, 0, 0), nrow = 8)
        b2 <- matrix(c(1, 0, 0, 0, 1, 1, 0, 0), nrow = 8)
        g1 <- 3
        g2 <- 3

    } else if (scenario == "7"){

        b1 <- matrix(c(1, 1, 0, 0, rep(0, 4)), nrow = 8)
        b2 <- matrix(c(1, 1, 0, 0, rep(0, 4)), nrow = 8)
        g1 <- 3
        g2 <- 3
    } else if (scenario == "8"){

        b1 <- matrix(c(1, 0, 0, .5, 0, 0, 0, .5), nrow = 8)
        b2 <- matrix(c(1, 0, 0, 0, rep(0, 4)), nrow = 8)
        g1 <- 3
        g2 <- 3
    } else stop("Scenario does not exist")


    list(b1 = b1, b2 = b2, g1 = g1, g2 = g2)

}


run_one_replicate <- function(scenario = "0"){

    b <- switch_scenario(scenario)

    data1out2 <- genset_func(n = 2500, beta1 = b$b1, gamma1 = b$g1, beta2 = b$b2, gamma2 = b$g2, scenario = scenario)

    dataou1 <- list(augment_data(data1out2[data1out2$trans == 1, ],
                                 data1out2$missind[data1out2$trans == 1], "ID", S ~ W, Nout=25),
                    augment_data(data1out2[data1out2$trans == 2, ],
                                 data1out2$missind[data1out2$trans == 2], "ID", S ~ W, Nout=25))

    fittedmods <- NULL
    for (i in 1:2) {
        fittedmods[[i]] <- tryCatch(
            estimate_flexparams(
                dataou1[[i]],
                "ID",
                Surv(time,  status) ~ W + V8 + S * tr,
                K = 4,
                method = "BFGS"
            ),
            error = function(e) {
                tryCatch(
                    estimate_flexparams(
                        dataou1[[i]],
                        "ID",
                        Surv(time,  status) ~ W + V8 + S * tr,
                        K = 3,
                        method = "BFGS"
                    ),
                    error = function(e) {
                        tryCatch(
                            estimate_flexparams(
                                dataou1[[i]],
                                "ID",
                                Surv(time,  status) ~ W + V8 + S * tr,
                                K = 2,
                                method = "BFGS"
                            ),
                            error = function(e) {
                                tryCatch(estimate_flexparams(
                                    dataou1[[i]],
                                    "ID",
                                    Surv(time,  status) ~ W + V8 + S * tr,
                                    K = 1,
                                    method = "BFGS"
                                ), error = function(e) {
                                    tryCatch(estimate_flexparams(
                                        dataou1[[i]],
                                        "ID",
                                        Surv(time,  status) ~ W + V8 + S * tr,
                                        K = 0,
                                        method = "BFGS"
                                    ), error = function(e) {
                                        list(optim = list(convergence = 99))
                                    })
                                })
                            }
                        )
                    }
                )

            }
        )
    }

    ## check convergence

    if(fittedmods[[1]]$optim$convergence == 0 & fittedmods[[2]]$optim$convergence == 0 ) {
        ##  post-estimation

        gcomp.est <- gcompute_cuminc(dataou1, fittedmods, timet = 2, "tr", "S", nn = 100)
        cuminc1.by.treatment <- tapply(gcomp.est$cuminc.cause1, list(gcomp.est$group, gcomp.est$tr), mean)
        cuminc2.by.treatment <- tapply(gcomp.est$cuminc.cause2, list(gcomp.est$group, gcomp.est$tr), mean)
        overall.surv.by.treatment <- tapply(gcomp.est$surv.overall, list(gcomp.est$group, gcomp.est$tr), mean)

        s.val <- tapply(gcomp.est$S, gcomp.est$group, mean)

        rdiff.1 <- (cuminc1.by.treatment[, 2] - cuminc1.by.treatment[, 1])
        rdiff.2 <- (cuminc2.by.treatment[, 2] - cuminc2.by.treatment[, 1])

        rd.s <- overall.surv.by.treatment[, 2] - overall.surv.by.treatment[, 1]


        list(convergence = c(fittedmods[[1]]$optim$convergence, fittedmods[[2]]$optim$convergence),
             s.vals = s.val,
             drd.s = rdiff.1 - rdiff.2,
             drd = mean(rdiff.1 - rdiff.2),
             rdiff.1 = rdiff.1,
             rdiff.2 = rdiff.2,
             rd.s = rd.s,
             fitted.objs = list(fittedmods[[1]], fittedmods[[2]])
        )

    } else {
        list(convergence = c(fittedmods[[1]]$optim$convergence, fittedmods[[2]]$optim$convergence),
             s.vals = NA,
             drd.s = NA,
             drd = NA,
             rdiff.1 = NA,
             rdiff.2 = NA,
             rd.s = NA,
             fitted.objs = list(fittedmods[[1]], fittedmods[[2]])
        )


    }

}
