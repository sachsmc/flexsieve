
set.seed(4402)
b <- switch_scenario("3")

data1out2 <- genset_func(n = 2500, beta1 = b$b1, gamma1 = b$g1, beta2 = b$b2, gamma2 = b$g2, scenario = "3")

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
                                K = 0,
                                method = "BFGS"
                            ), error = function(e) {
                                list(optim = list(convergence = 99))
                            })
                        }
                    )
                }
            )

        }
    )
}

saveRDS(fittedmods, file = "tests/example-fit.RData")


## bootstrap

bootmods <- vector(length = 200, mode = "list")

for(j in 1:200) {

bootid <- sample(unique(data1out2$ID), 2500, replace = TRUE)
data1out2b <- data1out2[c(sapply(bootid, function(x) which(data1out2$ID == x))), ]

dataou1b <- list(augment_data(data1out2b[data1out2b$trans == 1, ],
                             data1out2b$missind[data1out2b$trans == 1], "ID", S ~ W, Nout=25),
                augment_data(data1out2b[data1out2b$trans == 2, ],
                             data1out2b$missind[data1out2b$trans == 2], "ID", S ~ W, Nout=25))

fittedmods <- NULL
for (i in 1:2) {
    fittedmods[[i]] <- tryCatch(
        estimate_flexparams(
            dataou1b[[i]],
            "ID",
            Surv(time,  status) ~ W + V8 + S * tr,
            K = 4,
            method = "BFGS"
        ),
        error = function(e) {
            tryCatch(
                estimate_flexparams(
                    dataou1b[[i]],
                    "ID",
                    Surv(time,  status) ~ W + V8 + S * tr,
                    K = 3,
                    method = "BFGS"
                ),
                error = function(e) {
                    tryCatch(
                        estimate_flexparams(
                            dataou1b[[i]],
                            "ID",
                            Surv(time,  status) ~ W + V8 + S * tr,
                            K = 2,
                            method = "BFGS"
                        ),
                        error = function(e) {
                            tryCatch(estimate_flexparams(
                                dataou1b[[i]],
                                "ID",
                                Surv(time,  status) ~ W + V8 + S * tr,
                                K = 0,
                                method = "BFGS"
                            ), error = function(e) {
                                list(optim = list(convergence = 99))
                            })
                        }
                    )
                }
            )

        }
    )
}


bootmods[[j]] <- fittedmods

}


saveRDS(bootmods, file = "tests/bootstrap-fits.RData")


### analyze data

fittedmod <- readRDS("tests/example-fit.RData")
bootmods <- readRDS("tests/bootstrap-fits.RData")

## fix bootstraps that didn't converge

j.fail <- c(166, 184)

d1 <- NULL
d2 <- NULL
for(j in 1:200) {

bootid <- sample(unique(data1out2$ID), 2500, replace = TRUE)
data1out2b <- data1out2[c(sapply(bootid, function(x) which(data1out2$ID == x))), ]

dataou1b <- list(augment_data(data1out2b[data1out2b$trans == 1, ],
                             data1out2b$missind[data1out2b$trans == 1], "ID", S ~ W, Nout=25),
                augment_data(data1out2b[data1out2b$trans == 2, ],
                             data1out2b$missind[data1out2b$trans == 2], "ID", S ~ W, Nout=25))
if(j == j.fail[1]) d1 <- dataou1b
if(j == j.fail[2]) d2 <- dataou1b

}


bootmods.166.1 <- estimate_flexparams(
                                d1[[1]],
                                "ID",
                                Surv(time,  status) ~ W + V8 + S * tr,
                                K = 1,
                                method = "BFGS"
                            )

bootmods[[166]][[1]] <- bootmods.166.1



bootmods.184.1 <- estimate_flexparams(
                                d2[[1]],
                                "ID",
                                Surv(time,  status) ~ W + V8 + S * tr,
                                K = 2,
                                method = "BFGS"
                            )



bootmods[[j.fail[2]]][[1]] <- bootmods.184.1

## 4 panel figure

gcomp.est <- gcompute_cuminc(dataou1, fittedmod, timet = 2, "tr", "S", nn = 100)
cuminc1.by.treatment <- tapply(gcomp.est$cuminc.cause1, list(gcomp.est$group, gcomp.est$tr), mean)
cuminc2.by.treatment <- tapply(gcomp.est$cuminc.cause2, list(gcomp.est$group, gcomp.est$tr), mean)
overall.surv.by.treatment <- tapply(gcomp.est$surv.overall, list(gcomp.est$group, gcomp.est$tr), mean)

s.val <- tapply(gcomp.est$S, gcomp.est$group, mean)

rdiff.1 <- (cuminc1.by.treatment[, 2] - cuminc1.by.treatment[, 1])
rdiff.2 <- (cuminc2.by.treatment[, 2] - cuminc2.by.treatment[, 1])

rd.s <- -(overall.surv.by.treatment[, 2] - overall.surv.by.treatment[, 1])

mainres <- data.frame(bootdex = 0, s = s.val, drd.s = rdiff.1 - rdiff.2,
           rd.1 = rdiff.1, rd.2 = rdiff.2,
           rd.s = rd.s)


c1 <- fittedmod[[1]]$optim$par
c2 <- fittedmod[[2]]$optim$par

coeffs.main <- c(c1[(length(c1) - 1):length(c1)], c2[(length(c2) - 1):length(c2)])

coeffs.boot <- do.call(rbind, lapply(bootmods, function(fittedmod) {
    c1 <- fittedmod[[1]]$optim$par
    c2 <- fittedmod[[2]]$optim$par

    if(fittedmod[[1]]$optim$convergence != 0 |fittedmod[[2]]$optim$convergence != 0  ) return(rep(NA, 4))
    c(c1[(length(c1) - 1):length(c1)], c2[(length(c2) - 1):length(c2)])

}))

coeffs.boot <- coeffs.boot[!is.na(coeffs.boot[, 1]), ]
CIs <- apply(coeffs.boot, MAR = 2, FUN = function(x) sprintf("(%.2f to %.2f)", quantile(x, .025), quantile(x, .975)))
p.vals <- apply(coeffs.boot, MAR = 2, FUN = function(x) {
    2 * pnorm(-abs(median(x)) / (quantile(x, .75) - quantile(x, .25)))
})

xtable::xtable(cbind(round(coeffs.main, 2), CIs, scales::pvalue(p.vals)))

library(parallel)

cl <- makeCluster(getOption("cl.cores", 8))
clusterEvalQ(cl, devtools::load_all())
clusterExport(cl, varlist = c("bootmods"))
clusterExport(cl, varlist = c("dataou1"))

bootcurves <- clusterApplyLB(cl, 1:length(bootmods), function(i) {

    gcomp.est <- tryCatch(gcompute_cuminc(dataou1, bootmods[[i]], timet = 2, "tr", "S", nn = 50),
		error = function(e) NA)
	if(any(is.na(gcomp.est[[1]]))) {return(data.frame(bootdex = i, s = NA, drd.s = NA,
               rd.1 = NA, rd.2 = NA,
               rd.s = NA))
	}
    cuminc1.by.treatment <- tapply(gcomp.est$cuminc.cause1, list(gcomp.est$group, gcomp.est$tr), mean)
    cuminc2.by.treatment <- tapply(gcomp.est$cuminc.cause2, list(gcomp.est$group, gcomp.est$tr), mean)
    overall.surv.by.treatment <- tapply(gcomp.est$surv.overall, list(gcomp.est$group, gcomp.est$tr), mean)

    s.val <- tapply(gcomp.est$S, gcomp.est$group, mean)

    rdiff.1 <- (cuminc1.by.treatment[, 2] - cuminc1.by.treatment[, 1])
    rdiff.2 <- (cuminc2.by.treatment[, 2] - cuminc2.by.treatment[, 1])

    rd.s <- overall.surv.by.treatment[, 2] - overall.surv.by.treatment[, 1]

    data.frame(bootdex = i, s = s.val, drd.s = rdiff.1 - rdiff.2,
               rd.1 = rdiff.1, rd.2 = rdiff.2,
               rd.s = rd.s)

})

saveRDS(bootcurves, file = "tests/example-bootcurves.RData")

bootcurves <- readRDS("tests/example-bootcurves.RData")

library(data.table)
library(ggplot2)

bootplot <- data.table(do.call(rbind, bootcurves))
bootplot[, rd.s := -rd.s]

bootplot <- bootplot[!is.na(drd.s)]
bootquant <- bootplot[, .(lower.drd = quantile(drd.s, .025), upper.drd = quantile(drd.s, .975),
             lower.rd.s = quantile(rd.s, .025), upper.rd.s = quantile(rd.s, .975),
             lower.rd.1 = quantile(rd.1, .025), upper.rd.1 = quantile(rd.1, .975),
             lower.rd.2 = quantile(rd.2, .025), upper.rd.2 = quantile(rd.2, .975)), by = s]
toplot <- data.table(merge(bootquant, mainres, by = "s"))

p1 <- ggplot(toplot, aes(x = s, y = drd.s)) + geom_line() +
    geom_line(aes(y = lower.drd), linetype = 2) + scale_y_continuous(limits=c(-1, 1.5)) +
    geom_line(aes(y = upper.drd), linetype = 2) + theme_bw() + ggtitle("DRD curve") +
    xlab(expression(s[1])) + ylab(expression(DRD(s[1]))) + geom_rug(sides = "b")

p2 <- ggplot(toplot, aes(x = s, y = rd.s)) + geom_line() +
    geom_line(aes(y = lower.rd.s), linetype = 2) + scale_y_continuous(limits=c(-1, 1.5)) +
    geom_line(aes(y = upper.rd.s), linetype = 2)+ theme_bw() + ggtitle("RD curve") +
    xlab(expression(s[1])) + ylab(expression(RD(s[1]))) + geom_rug(sides = "b")

p3 <- ggplot(toplot, aes(x = s, y = rd.1)) + geom_line() +
    geom_line(aes(y = lower.rd.1), linetype = 2) + scale_y_continuous(limits=c(-1, 1.5)) +
    geom_line(aes(y = upper.rd.1), linetype = 2)+ theme_bw() + ggtitle("RD(s) matched")

p4 <- ggplot(toplot, aes(x = s, y = rd.2)) + geom_line() +
    geom_line(aes(y = lower.rd.2), linetype = 2) + scale_y_continuous(limits=c(-1, 1.5)) +
    geom_line(aes(y = upper.rd.2), linetype = 2)+ theme_bw() + ggtitle("RD(s) mismatched")

library(patchwork)

p1 + p2
ggsave("example-figure.pdf", width = 6, height = 3.5)


bootplot[, .(m.drd = mean(drd.s)), by = bootdex][, .(lower = quantile(m.drd, .025), upper= quantile(m.drd, .975))]
mean(mainres$drd.s)
