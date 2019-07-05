
#' Augment data frame for missing counterfactuals
#'
#' Impute missing counterfactuals using a linear model on the
#' baseline measurements
#'
#' @param data1 Data frame in question
#' @param missingind vector the same number of rows as data1, indicator for the missingness
#' @param ID.var character vector specifying the ID variable in data1
#' @param formula.predict Formula specifying the model
#' @param Nout integer, the number of imputed observations per missing one
#'
#' @return A list containing two dataframes, the unaugmented data, and the augmented data
#'
#' @export

augment_data<-function(data1, missingind, ID.var, formula.predict, Nout){

    datavad<-data1[missingind!=1,]
    datamid<-data1[missingind==1,]

    predfit <- lm(formula.predict, data = datavad)
    modout<-summary(predfit)

    datamid<-datamid[order(datamid[[ID.var]]),]
    SQhatmat<- c(model.matrix(update(formula.predict, NULL ~ .), data = datamid) %*% modout$coefficients[,1, drop = FALSE])
    s1mat<-unlist(lapply(SQhatmat,function(x){rnorm(Nout,x,modout$sigma)}))
    bigdataout<-datamid
    for(i in 1:(Nout-1)){
        bigdataout<-rbind(bigdataout,datamid)
    }
    bigdataout<-bigdataout[order(bigdataout[[ID.var]]),]

    key1 <- all.vars(formula.predict)[1]
    bigdataout[[key1]]<-s1mat


    return(list(data1=datavad, data2=bigdataout,
                fitted.w.s = predfit))
}


#' Estimate parameters in the flexsurv model
#'
#' @param datain List of data frames as returned by \link{augment_data}
#' @param ID.var characted indicating id variable in data
#' @param form Formula specifying the model. The outcome should be a \link{Surv} object
#' @param K the number of knots
#' @param ... Other arguments passed to \link{optim}
#'
#' @export

estimate_flexparams <- function(datain, ID.var, form ,K, ...){

    data2<-datain$data2
    data1<-datain$data1
    IDout <- data2[[ID.var]]

    lhood2 <- flexsurvspline.likelihood(form, data = data2, k = K) ## function of parameters
    lhood1 <- flexsurvspline.likelihood(form, data = data1, k = K) ## function of parameters

    nparm <- ncol(model.matrix(form, data1)) - 1 + (K + 2)

    sfit <- tryCatch(flexsurvspline(form, data = rbind(data1, data2), k = K),
                     error = function(e) list(coefficients = rep(0.01, nparm)))
    start <- sfit$coefficients

    lilwyan<-function(beta1){

        out1<-lhood2(beta1)
        out2<-pmax(1e-12, tapply(out1, IDout, mean))

        loglikout<-sum(log(pmax(1e-12, lhood1(beta1))))+sum(log(out2))

        return(-loglikout)
    }


    beta0<- start

    optimout1<-optim(beta0,lilwyan, ...)

    sfit$res[, 1] <- sfit$coefficients <- optimout1$par

    retlist <- list(optim = optimout1,
                    knots = sfit$knots,
                    formula = update(form, ~ . -1),
                    flexobj = sfit)
    return(retlist)

}

