library(flexsieve)

## setup test data
data(lung)

data1out<- lung[sample(1:nrow(lung), 400, replace = TRUE), ]
data1out$treatment <- rbinom(dim(data1out)[1],1,0.3)
data1out$missind<- 1 - data1out$treatment
data1out$ID<-1:length(data1out$missind)
data1out$pat.karno <- rnorm(400)
data1out$ph.karno <- rnorm(400)
data1out$statusm <- data1out$status - 1
data1out$statusm[sample(1:nrow(data1out))[1:80]] <- 2
data1out$status1 <- as.numeric(data1out$statusm == 1)
data1out$status2 <- as.numeric(data1out$statusm == 2)
data1out$start <- 0
data1out$status0 <- 1

Tmat <- rbind(c(NA, 1, 2), c(NA, NA, NA), c(NA, NA, NA))

data1out2 <- msprep(c("start", "time", "time"), c("status0", "status1", "status2"), data1out, trans = Tmat,
                    id = "ID", keep = c("age", "pat.karno", "ph.karno", "treatment", "missind"))
## main functions in the package

dataou1 <- list(augment_data(data1out2[data1out2$trans == 1, ],data1out2$missind[data1out2$trans == 1], "ID", ph.karno ~ pat.karno, Nout=25),
                augment_data(data1out2[data1out2$trans == 2, ],data1out2$missind[data1out2$trans == 2], "ID", ph.karno ~ pat.karno, Nout=25))

fittedmods <- NULL
for(i in 1:2) {
    fittedmods[[i]] <- estimate_flexparams(dataou1[[i]], "ID", Surv(time,  status) ~ age + pat.karno + ph.karno * treatment, K = 3, method = "BFGS")
}

system.time(gcomp.est <- gcompute_cuminc(dataou1, fittedmods, timet = 450, "treatment", "ph.karno", nn = 40))

## summarize the results

cuminc.by.treatment <- tapply(gcomp.est$cuminc, list(gcomp.est$group, gcomp.est$treatment), mean)
s.val <- tapply(gcomp.est$ph.karno, gcomp.est$group, mean)

plot(cuminc.by.treatment[, 1] ~ s.val, ylim = c(0, 1))
lines(cuminc.by.treatment[, 2] ~ s.val, col = "red", type = "p")
