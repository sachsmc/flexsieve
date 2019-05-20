
data(lung)

data1out<-lung
data1out$treatment <- rbinom(dim(lung)[1],1,0.3)
data1out$missind<- data1out$treatment
data1out$ID<-1:length(data1out$missind)
data1out$pat.karno[is.na(data1out$pat.karno)]<-1
data1out$ph.karno[is.na(data1out$ph.karno)]<-1

dataou1 <- augment_data(data1out,data1out$missind, "ID", ph.karno ~ pat.karno, Nout=25)
fittedmod <- estimate_flexparams(dataou1, "ID", Surv(time,  status == 1) ~ age + pat.karno + ph.karno * treatment, K = 3, method = "BFGS")

gcomp.est <- gcompute_cuminc(dataou1, fittedmod, timet = 450, "treatment", "ph.karno", nn = 100)

cuminc.by.treatment <- tapply(gcomp.est$cuminc, list(gcomp.est$group, gcomp.est$treatment), mean)
s.val <- tapply(gcomp.est$ph.karno, gcomp.est$group, mean)

plot(cuminc.by.treatment[, 1] ~ s.val, ylim = c(0, 1))
lines(cuminc.by.treatment[, 2] ~ s.val, col = "red", type = "p")
