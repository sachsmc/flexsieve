# cbind(1,tr,S1,tr*S1,W,X3)
check_conv <- function(i) {


    b1 <- matrix(c(1.3, 0, 0, .5, rep(0, 4)), nrow = 8)
    b2 <- matrix(c(1.3, 0, 0, -.5, rep(0, 4)), nrow = 8)
    g1 <- 3
    g2 <- 3

data1out2 <- genset_func(n = 2000, beta1 = b1 ,gamma1 = g1, beta2 = b2, gamma2 = g2)

dataou1 <- list(augment_data(data1out2[data1out2$trans == 1, ],
                             data1out2$missind[data1out2$trans == 1], "ID", S ~ W, Nout=25),
                augment_data(data1out2[data1out2$trans == 2, ],
                             data1out2$missind[data1out2$trans == 2], "ID", S ~ W, Nout=25))

fittedmods <- NULL
formu <- Surv(time,  status) ~ tr + S + tr:S + V8 + W
for(i in 1:2) {
    fittedmods[[i]] <- tryCatch(estimate_flexparams(dataou1[[i]], "ID",
                                           formu,
                                           K = 4, method = "BFGS"),
                                error = function(e) {
                                    tryCatch(estimate_flexparams(dataou1[[i]], "ID",
                                                                 formu,
                                                                 K = 3, method = "BFGS"), error = function(e) {
                                                                     tryCatch(estimate_flexparams(dataou1[[i]], "ID",
                                                                                         formu,
                                                                                         K = 2, method = "BFGS"),
                                                                              error = function(e) {
                                                                                  estimate_flexparams(dataou1[[i]], "ID",
                                                                                                      formu,
                                                                                                      K = 0, method = "BFGS")
                                                                              })
                                                                 })

                                })
}

c(fittedmods[[1]]$optim$convergence, fittedmods[[1]]$optim$par, fittedmods[[2]]$optim$convergence, fittedmods[[2]]$optim$par)

}

lillone <- sapply(1:10, check_conv)

res2 <- do.call(rbind, lapply(lillone, function(x){
    x[!grepl("gamma", names(x))]
}))

#gcomp.est <- gcompute_cuminc(dataou1, fittedmods, timet = 2, "tr", "S", nn = 40)

cuminc.by.treatment <- tapply(gcomp.est$cuminc.cause1, list(gcomp.est$group, gcomp.est$tr), mean)
s.val <- tapply(gcomp.est$S, gcomp.est$group, mean)

plot(cuminc.by.treatment[, 1] ~ s.val, ylim = c(0, 1))
lines(cuminc.by.treatment[, 2] ~ s.val, col = "red", type = "p")



system.time(test1 <- run_one_replicate("0"))
