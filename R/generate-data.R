#' Generate survival data for analysis
#'
#' 2500 subjects, 20 variables, about XX% cumulative incidence. 5 variables truly associated,
#' some linear, some wonky associations, an interaction, and some variables correlated.
#'
#' n Sample size
#' Character designating the simulation scenario:
#'                 0 = NULL, no association,
#'                 A = simple, one variable is linearly associated with IRIS,
#'                 B = spline of 5
#'                 C = wonky, interactions, nonlinearities, etc.
#'                 D = same as C
#' Proportion of missing binary Y values (random censoring)
#'
#'A data frame with X variables, censored survival times (competing risk w death), and true cumulative incidence at Tstar weeks
#'
#'

Haz11<-function(t,Z, beta,gamma){

    lambda <- exp(Z %*% beta)
    (t / lambda)^gamma

}

haz11 <- function(t, Z, beta, gamma) {

    lambda <- exp(Z %*% beta)
    gamma * t^(gamma - 1) * (1/lambda)^gamma

}



St <- function(t, Z,beta1,gamma1,beta2,gamma2){

exp(-(Haz11(t,Z, beta1,gamma1)+Haz11(t,Z, beta2,gamma2)))

}



inverse <- function(fn, min_x, max_x,p) {

    fn_inv = function(y){
        uniroot((function(x){fn(x) - y}), lower=min_x, upper=max_x)[1]$root
    }
    return(Vectorize(fn_inv))
}


sampletimesfunc <- function(mat,b1,g1,b2,g2) {
    n <- dim(mat)[1]
    UU <- runif(n)

    sapply(1:n, function(i){

        inverse(function(t) St(t, mat[i,], beta1 = b1,gamma1=g1, beta2=b2, gamma2=g2), 0, 1e8)(UU[i])

    })
}





genset_func<-function(n, beta1, gamma1, beta2, gamma2, scenario){

    tr<-rbinom(n,1,0.5)
    X3 <- matrix(rnorm(n * 3, mean = 0, sd = .1), ncol = 3)
    W<-rnorm(n,2,0.5)

    if(scenario == "4") {
        S1 <- rnorm(n, 0.3 + .33 * W, 0.75)
    } else {
        S1<-rnorm(n,0.3+W,0.75)
    }
    S1 <- S1 - mean(S1)

    if(scenario == "6") { # dependent censoring

        Cen <- rexp(n, rate= 1 / exp(X3[, 1] * .5 + X3[, 2] * .5))

    } else if (scenario == "7") {

        Cen <- rexp(n, rate = 1 / exp(S1 * .5))

    } else {

        p2 <- runif(n)
        Cen<- exp(log(-log(1-p2))+1.5)

    }

    if(!is.matrix(beta1)) {
    beta100<-c(beta1, rep(0,8-length(beta1)))
    beta10<-as.matrix(beta100,nrow=8)

    beta200<-c(beta2, rep(0,8-length(beta2)))
    beta20<-as.matrix(beta100,nrow=8)
    } else {
        beta10 <- beta1
        beta20 <- beta2
    }
    matx<-cbind(1,tr,S1,tr*S1,X3, W)
    Touts<-sampletimesfunc(matx,beta10,gamma1,beta20,gamma2)


        type1<-rbinom(n,1,haz11(Touts,matx,beta10,gamma1)/(haz11(Touts,matx,beta10,gamma1)+haz11(Touts,matx,beta20,gamma2)))



    type1<-ifelse(type1 == 1, 1, 2)

    ID<-1:length(type1)

    Toutfin<-pmin(Touts, Cen)

    type1<-ifelse(Cen<Touts, 0, type1)

    S<-ifelse(tr==1, S1, NA)
    dat1<-as.data.frame(cbind(ID, Toutfin, type1, tr, S,S1, W, X3))

    dat1$status1 <- as.numeric(dat1$type1 == 1)
    dat1$status2 <- as.numeric(dat1$type1 == 2)
    dat1$start <- 0
    dat1$status0 <- 1
    dat1$missind <- as.numeric(is.na(dat1$S))

    Tmat <- rbind(c(NA, 1, 2), c(NA, NA, NA), c(NA, NA, NA))

    data1out2 <- msprep(c("start", "Toutfin", "Toutfin"), c("status0", "status1", "status2"), dat1, trans = Tmat,
                        id = "ID", keep = c("missind", "tr", "S", "W", "V8", "V9", "V10"))
    ## main functions in the package

    return(data1out2)
}
