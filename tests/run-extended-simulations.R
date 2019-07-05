library(parallel)

cl <- makeCluster(getOption("cl.cores", 8))
clusterEvalQ(cl, devtools::load_all())

ree.scen <- function(scen) {

    res.X <- clusterApplyLB(cl, 1:200, function(x){
        res.1 <- run_one_replicate(scen)
        saveRDS(res.1, file =
                    file.path("dataraw",
                              sprintf("sim-res-%s-%s.RData", scen,
                                      formatC(x, width = 3, format = "d", flag = "0"))))
    })

}

ree.scen("4")
ree.scen("5")
ree.scen("6")
ree.scen("7")
ree.scen("8")


pull.scen <- function(scenario) {

    files <- unlist(lapply(1:200, function(x) {
        fin <- file.path("dataraw",
                         sprintf("sim-res-%s-%s.RData", scenario,
                                 formatC(x, width = 3, format = "d", flag = "0")))
        if(file.exists(fin)) return(fin) else return(NULL)
    }))

    b <- switch_scenario(scenario)

    truestuff <- true.curves(b$b1, b$b2, b$g1, b$g2, scenario)

    simres <- lapply(1:length(files), function(i) {
        res0 <- readRDS(files[i])
        data.frame(scenario = scenario, replicate = i, converged = all(res0$convergence == 0),
                   s = res0$s.vals, drd.s = res0$drd.s, rd.1 = res0$rdiff.1, rd.2 = res0$rdiff.2,
                   rd.s = -res0$rd.s)
    })
    simout <- do.call(rbind, simres)
    list(simout, truestuff)

}


scen4 <- pull.scen("4")
scen5 <- pull.scen("5")
scen6 <- pull.scen("6")
scen7 <- pull.scen("7")
scen8 <- pull.scen("8")

truescens <- rbind(data.frame(scenario = "0", replicate = 0, converged = TRUE,
                              s = scen0[[2]]$s, drd.s = scen0[[2]]$drd.s, rd.s = scen0[[2]]$rd.curve),
                   data.frame(scenario = "A", replicate = 0, converged = TRUE,
                              s = scen1[[2]]$s, drd.s = scen1[[2]]$drd.s, rd.s = scen1[[2]]$rd.curve),
                   data.frame(scenario = "B", replicate = 0, converged = TRUE,
                              s = scen2[[2]]$s, drd.s = scen2[[2]]$drd.s, rd.s = scen2[[2]]$rd.curve),
                   data.frame(scenario = "C", replicate = 0, converged = TRUE,
                              s = scen3[[2]]$s, drd.s = scen3[[2]]$drd.s, rd.s = scen3[[2]]$rd.curve))


library(ggplot2)
library(data.table)
library(quantreg)
library(mgcv)

msd <- function(x, ...) {
    q3 <- quantile(x, .75, ...)
    q1 <- quantile(x, .25, ...)
    abs(q3 - q1)
}


allscens <- as.data.table(rbind(scen0[[1]], scen1[[1]], scen2[[1]], scen3[[1]]))
allscens[, scenario := c("0", "A", "B", "C")[as.numeric(scenario)]]

allscens[, s.cat := cut(s, breaks = seq(-4, 4, by = .1), labels = FALSE)]
allscens[, s.cat.val := (seq(-4, 4, by = .1) - .05)[s.cat]]

qnat <- allscens[, .(upper.95.drd = quantile(drd.s, .975, na.rm = TRUE),
                     lower.95.drd = quantile(drd.s, .025, na.rm = TRUE),
                     upper.95.rd = quantile(rd.s, .975, na.rm = TRUE),
                     lower.95.rd = quantile(rd.s, .025, na.rm = TRUE),
                     sd.drd = sd(drd.s, na.rm = TRUE),
                     sd.rd = sd(rd.s, na.rm = TRUE)), by = .(scenario, s.cat.val)]

allscens2 <- merge(allscens,  qnat[, .(scenario, s.cat.val, sd.drd, sd.rd)], by = c("scenario", "s.cat.val"))
allscens2 <- allscens2[!is.na(s.cat) & converged == TRUE]

ggplot(allscens, aes(x = s, y = drd.s)) +
    geom_smooth(se = FALSE) +
    geom_line(data = qnat, aes(x = s.cat.val, y = upper.95.drd), color = "blue", linetype = 3) +
    geom_line(data = qnat, aes(x = s.cat.val, y = lower.95.drd), color = "blue", linetype = 3) +
    geom_line(data = truescens, aes(x = s, y = drd.s), color = "red") +
    facet_wrap(~ scenario) + ggtitle("DRD curve") + theme_bw()
ggsave(file ="drd-curves.png", width = 6.5, height = 5.5)


ggplot(allscens, aes(x = s, y = rd.s)) +
    geom_smooth(se = FALSE) +
    geom_line(data = qnat, aes(x = s.cat.val, y = upper.95.rd), color = "blue", linetype = 3) +
    geom_line(data = qnat, aes(x = s.cat.val, y = lower.95.rd), color = "blue", linetype = 3) +
    geom_line(data = truescens, aes(x = s, y = rd.s), color = "red") +
    facet_wrap(~ scenario) + ggtitle("Risk difference") + theme_bw() + ggtitle("RD curve")
ggsave(file ="rd-curves.png", width = 6.5, height = 5.5)



marg.drd <- allscens[, .(mean.drd = mean(drd.s)), by = .(scenario, replicate)][, .(mean.drd = mean(mean.drd, na.rm = TRUE),
                                                                                   sd.drd = sd(mean.drd, na.rm = TRUE)), by = scenario]

marg.drd[, true.drd := c(scen0[[2]]$drd, scen1[[2]]$drd, scen2[[2]]$drd, scen3[[2]]$drd)]


## testable hypotheses

truestuff <- as.data.table(truescens)
truestuff[, ":="(s.cat.val = s, true.drd = drd.s, true.rd = rd.s)]

allscens3 <- merge(allscens2, truestuff[, .(scenario, s.cat.val, true.drd, true.rd)], by = c("scenario", "s.cat.val"))

## drd excludes 0 at any point
# compute new critical value
library(mvtnorm)

q.band <- rep(NA, 4)
names(q.band) <- c("0", "A", "B", "C")

for(scen in c("0", "A", "B", "C")) {
    dmat <- cov(as.matrix(dcast(allscens2[scenario == scen & abs(s < 2), .(drd.s), keyby = .(replicate, s.cat.val)],
                                replicate ~ s.cat.val, value.var = "drd.s", fun.aggregate= mean, na.rm = TRUE)[, -1]), use = "pairwise")
    dmat[is.na(dmat)] <- 0

    samp.1 <- rmvnorm(500, mean = rep(0, nrow(dmat)), sigma = dmat)
    q.band[scen] <- quantile(apply(abs(samp.1) / t(sapply(1:500, function(i) pmax(sqrt(diag(dmat)), .01))), MAR = 1, FUN = max), .95)
}

allscens2[, q.band := q.band[scenario]]

allscens2[, reject.drd := sign(drd.s - q.band * sd.drd) == sign(drd.s + q.band * sd.drd)]
allscens2[, reject.rd := sign(rd.s - q.band * sd.rd) == sign(rd.s + q.band * sd.rd)]

pow.drd.s <- allscens2[, .(any.s.rej = any(reject.drd),
                           any.s.rd = any(reject.rd)),
                       by = .(scenario, replicate)][, .(power.drd.s = mean(any.s.rej),
                                                        power.rd.s = mean(any.s.rd)), keyby = scenario]


## total gain for drd

est_tg <- function(s, y) {

    keep <- !duplicated(s)
    s1 <- sort(s[keep])
    y1 <- y[keep][order(s[keep])]
    n <- length(s1)
    widths <- s1[2:n] - s1[1:(n - 1)]
    heights <- abs((y1[2:n] + y1[1:(n - 1)]) / 2 - mean(y, na.rm = TRUE))

    sum(widths * heights)

}

total_gains <- allscens[converged == TRUE, .(tg_drd = est_tg(s, drd.s), tg_rd = est_tg(s, rd.s),
                                             tg_rd1 = est_tg(s, rd.1), tg_rd2 = est_tg(s, rd.2)), by = .(scenario, replicate)]

total_gains[, ":="(sd.tg.drd = sd(tg_drd), sd.tg.rd = sd(tg_rd),
                   sd.tg.rd1 = sd(tg_rd1), sd.tg.rd2 = sd(tg_rd2)), by = scenario]

total_gains[, ":="(reject.tg.drd = sign(tg_drd - 1.96 * sd.tg.drd) == sign(tg_drd + 1.96 * sd.tg.drd),
                   reject.tg.rd = sign(tg_rd - 1.96 * sd.tg.rd) == sign(tg_rd + 1.96 * sd.tg.rd),
                   reject.tg.rd1 = sign(tg_rd1 - 1.96 * sd.tg.rd1) == sign(tg_rd1 + 1.96 * sd.tg.rd1),
                   reject.tg.rd2 = sign(tg_rd2 - 1.96 * sd.tg.rd2) == sign(tg_rd2 + 1.96 * sd.tg.rd2))]

tg_power <- total_gains[, .(pow.tg.drd = mean(reject.tg.drd),
                            pow.tg.rd = mean(reject.tg.rd),
                            pow.tg.rd1 = mean(reject.tg.rd1),
                            pow.tg.rd2 = mean(reject.tg.rd2)), by = scenario]

### slope based testable hypothesis

est_slope <- function(s, y) {

    rise <- y[which.max(s)] - y[which.min(s)]
    run <- max(s) - min(s)

    rise / run

}

slopes <- allscens[converged == TRUE, .(sl_drd = est_slope(s, drd.s), sl_rd = est_slope(s, rd.s),
                                        sl_rd1 = est_slope(s, rd.1), sl_rd2 = est_slope(s, rd.2)), by = .(scenario, replicate)]

slopes[, ":="(sd.sl.drd = sd(sl_drd), sd.sl.rd = sd(sl_rd), sd.sl.rd1 = sd(sl_rd1), sd.sl.rd2 = sd(sl_rd2)), by = scenario]
slopes[, ":="(rej.sl.drd = sign(sl_drd - 1.96 * sd.sl.drd) == sign(sl_drd + 1.96 * sd.sl.drd),
              rej.sl.rd = sign(sl_rd - 1.96 * sd.sl.rd) == sign(sl_rd + 1.96 * sd.sl.rd),
              rej.sl.rd1 = sign(sl_rd1 - 1.96 * sd.sl.rd1) == sign(sl_rd1 + 1.96 * sd.sl.rd1),
              rej.sl.rd2 = sign(sl_rd2 - 1.96 * sd.sl.rd2) == sign(sl_rd2 + 1.96 * sd.sl.rd2))]

slope.power <- slopes[, .(pow.sl.drd = mean(rej.sl.drd), pow.sl.rd = mean(rej.sl.rd),
                          pow.sl.rd1 = mean(rej.sl.rd1), pow.sl.rd2 = mean(rej.sl.rd2)), by = scenario]

round(t(slope.power[, -1]), 2)

##  marginal drd

marg.drd.test <- allscens[, .(marg.drd = mean(drd.s), marg.rd = mean(rd.s)), by = .(scenario, replicate)][!is.na(marg.drd)]
marg.drd.test[, sd.drd := msd(marg.drd), by = scenario]
marg.drd.test[, sd.rd := msd(marg.rd), by = scenario]

marg.drd.test[, reject.drd := sign(marg.drd -1.96 * sd.drd) == sign(marg.drd + 1.96 * sd.drd)]
marg.drd.test[, reject.rd := sign(marg.rd -1.96 * sd.rd) == sign(marg.rd + 1.96 * sd.rd)]

pow.drd.marg <- marg.drd.test[, .(power.drd = mean(reject.drd),
                                  power.rd = mean(reject.rd)), by = scenario]


####


pull.coeff <- function(scenario) {

    files <- unlist(lapply(1:200, function(x) {
        fin <- file.path("dataraw",
                         sprintf("sim-res-%s-%s.RData", scenario,
                                 formatC(x, width = 3, format = "d", flag = "0")))
        if(file.exists(fin)) return(fin) else return(NULL)
    }))




    simres <- lapply(1:length(files), function(i) {
        res0 <- readRDS(files[i])
        conv <- res0$convergence
        if(!all(conv == 0)) return(rep(NA, 13))
        par1 <- res0$fitted.objs[[1]]$optim$par
        par2 <- res0$fitted.objs[[2]]$optim$par
        c(replicate = i, conv[1], par1[(length(par1) - 4):length(par1)], conv[2], par2[(length(par2) - 4):length(par2)])

    })
    simout <- as.data.frame(do.call(rbind, simres))
    simout$scenario <- scenario
    simout

}

##  coefficient for Z in each cause
##  coefficient for interaction term

allcoeffs <- data.table(do.call(rbind, lapply(c("0", "1", "2", "3"), pull.coeff)))
colnames(allcoeffs) <- c("replicate", paste(c("converge", "W", "V8", "S", "tr", "S.tr"), rep(c("c1", "c2"), each = 6), sep = "."), "scenario")

plot(rd.s ~ s, data = allscens[scenario == "A" & replicate == 163])
lines(-rd.s ~ s, data = subset(truescens, scenario == "A"), col = "red")

ggplot(allscens[scenario == "A" & replicate == 163], aes(x = s, y = rd.s)) + geom_line() +
    geom_line(data = allscens3[scenario == "A" & replicate == 1], aes(x = s.cat.val, y = true.rd), color = "red")

allcoeffs <- allcoeffs[!is.na(converge.c1)]


allcoeffs[, ":="(sd.tr.c1 = msd(tr.c1), sd.tr.c2 = msd(tr.c2),
                 sd.inter.c1 = msd(S.tr.c1), sd.inter.c2 = msd(S.tr.c2)), by = scenario]

testme <- function(b, sd.b) {

    sign(b -1.96 * sd.b) == sign(b + 1.96 * sd.b)

}

allcoeffs[, ":="(test.z1 = testme(tr.c1, sd.tr.c1),
                 test.z2 = testme(tr.c2, sd.tr.c2),
                 test.int1 = testme(S.tr.c1, sd.inter.c1),
                 test.int2 = testme(S.tr.c2, sd.inter.c2))]

pow.coeff <- allcoeffs[, .(pow.z1 = mean(test.z1), pow.z2 = mean(test.z2),
                           pow.inter1 = mean(test.int1), pow.inter2 = mean(test.int2)), by = scenario]

pow.coeff[, scenario := c("0", "A", "B", "C")[as.numeric(scenario) + 1]]

xtable::xtable(t(allcoeffs[, .(sprintf("%.2f (%.2f)", median(tr.c1), sd.tr.c1[1]),
                               sprintf("%.2f (%.2f)", median(tr.c2), sd.tr.c2[1]),
                               sprintf("%.2f (%.2f)", median(S.tr.c1), sd.inter.c1[1]),
                               sprintf("%.2f (%.2f)", median(S.tr.c2), sd.inter.c2[1])), by = scenario]))


power.table <- merge(merge(pow.drd.s, pow.drd.marg, by = "scenario"), pow.coeff, by = "scenario")

pt2 <- t(power.table[, -1])
colnames(pt2) <- c("0", "A", "B", "C")
xtable::xtable(pt2)


mean(allcoeffs[scenario == "1"]$S.tr.c1)
mean(allcoeffs[scenario == "1"]$S.tr.c2)
