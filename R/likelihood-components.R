logLikFactory2 <- function(Y, X=0, weights, bhazard, dlist,inits, dfns, aux, mx, fixedpars=NULL) {
    pars   <- inits
    npars  <- length(pars)
    #message("N parameters: ", npars, "\n")
    nbpars <- length(dlist$pars)
    insert.locations <- setdiff(seq_len(npars),
                                fixedpars)

    ## which are the subjects with events
    event <- Y[,"status"] == 1
    event.times <- Y[event, "time1"]
    left.censor <- Y[!event, "time2"]
    right.censor <- Y[!event, "time1"]

    event.weights <- weights[event]
    no.event.weights <- weights[!event]

    par.transform <- flexsurv:::buildTransformer(inits, nbpars, dlist)

    aux.pars <- flexsurv:::buildAuxParms(aux, dlist)

    default.offset <- rep.int(1, length(event.times))
    do.hazard <- any(bhazard > 0)

    lik <- rep.int(0, nrow(Y))
    ## the ... here is to work around optim
    function(optpars, ...) {

        pars[insert.locations] <- optpars
        raw.pars <- pars
        pars <- as.list(pars)

        pars.event <- pars.nevent <- pars
        if (npars > nbpars) {
            beta <- raw.pars[(nbpars+1):npars]
            for (i in dlist$pars){
                pars[[i]] <- pars[[i]] +
                    X[,mx[[i]],drop=FALSE] %*% beta[mx[[i]]]
                pars.event[[i]] <- pars[[i]][event]
                pars.nevent[[i]] <- pars[[i]][!event]
            }
        }

        fnargs <- c(par.transform(pars),
                    aux.pars)
        fnargs.event <- c(par.transform(pars.event),
                          aux.pars)
        fnargs.nevent <- c(par.transform(pars.nevent),
                           aux.pars)

        ## Generic survival model likelihood contributions
        ## Observed deaths
        dargs <- fnargs.event
        dargs$x <- event.times
        dargs$log <- FALSE
        dens <- do.call(dfns$d, dargs)

        ## Left censoring times (upper bound for event time)
        if (any(!event)){
            pmaxargs <- fnargs.nevent
            pmaxargs$q <- left.censor # Inf if right-censored, giving pmax=1
            pmax <- do.call(dfns$p, pmaxargs)
            pmax[pmaxargs$q==Inf] <- 1  # in case user-defined function doesn't already do this

            ## Right censoring times (lower bound for event time)
            pargs <- fnargs.nevent
            pargs$q <- right.censor
            pmin <-  do.call(dfns$p, pargs)
        }

        ## Left-truncation
        targs   <- fnargs
        targs$q <- Y[,"start"]
        pobs <- 1 - do.call(dfns$p, targs) # prob of being observed = 1 unless left-truncated

        ## Hazard offset for relative survival models
        if (do.hazard){
            pargs   <- fnargs.event
            pargs$q <- event.times
            pminb   <- do.call(dfns$p, pargs)
            haz  <- dens / (1 - pminb)
            offseti <- (1 + bhazard[event] / (haz)*weights[event])
        } else {
            offseti <- default.offset
        }
        ## Express as vector of individual likelihood contributions
        lik[event] <- (dens*event.weights)
        if (any(!event))
            lik[!event] <- (pmax - pmin)*no.event.weights
        lik / pobs
    }
}


flexsurvreg.likelihood <- function(formula, anc=NULL, data, weights, bhazard, subset, na.action, dist,
                                   inits, fixedpars=NULL, dfns=NULL, aux=NULL, cl=0.95,
                                   integ.opts=NULL, sr.control=survreg.control(), ...)
{
    call <- match.call()
    if (missing(dist)) stop("Distribution \"dist\" not specified")
    if (is.character(dist)) {
        # Added: case insensitve matching of distributions
        # Step 1: Use match.arg on lowercase argument, dists.
        # Step 2: Use match to get index in distribution list from
        # value obtained in step 1, and grab corresponding element.
        dist <- match.arg(tolower(dist), tolower(names(flexsurv.dists)))
        dist <- names(flexsurv.dists)[match(dist,tolower(names(flexsurv.dists)))]
        dlist <- flexsurv.dists[[dist]]
    }
    else if (is.list(dist)) {
        dlist <- flexsurv:::check.dlist(dist)
    }
    else stop("\"dist\" should be a string for a built-in distribution, or a list for a custom distribution")
    dfns <- flexsurv:::form.dp(dlist, dfns, integ.opts)
    parnames <- dlist$pars
    ancnames <- setdiff(parnames, dlist$location)

    flexsurv:::check.formula(formula, dlist)
    if (is.null(anc)){
        anc <- vector(mode="list", length=length(ancnames))
        names(anc) <- ancnames
        for (i in ancnames){
            anc[[i]] <- flexsurv:::ancpar.formula(formula, i)
        }
    }
    else {
        if (!is.list(anc) || !all(sapply(anc, function(x)inherits(x, "formula"))))
            stop("\"anc\" must be a list of formulae")
    }
    forms <- c(location=flexsurv:::get.locform(formula, ancnames), anc)
    names(forms)[[1]] <- dlist$location

    ## a) calling model.frame() directly doesn't work.  it only looks in
    ## "data" or the environment of "formula" for extra variables like
    ## "weights". needs to also look in environment of flexsurvreg.
    ## m <- model.frame(formula=, data=data, weights=weights, subset=subset, na.action=na.action)
    ## b) putting block below in a function doesn't work when calling
    ## flexsurvreg within a function
    ## m <- make.model.frame(call, formula, data, weights, subset, na.action, ancnames)

    ## Make model frame
    indx <- match(c("formula", "data", "weights", "bhazard", "subset", "na.action"), names(call), nomatch = 0)
    if (indx[1] == 0)
        stop("A \"formula\" argument is required")
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")

    f2 <- flexsurv:::concat.formulae(formula,forms)
    temp[["formula"]] <- f2
    if (missing(data)) temp[["data"]] <- environment(formula)
    m <- eval(temp, parent.frame())
    m <- droplevels(m) # remove unused factor levels after subset applied
    attr(m,"covnames") <- attr(f2, "covnames") # for "newdata" in summary
    attr(m,"covnames.orig") <- intersect(colnames(m), attr(f2, "covnames.orig")) # for finding factors in plot method
    Y <- flexsurv:::check.flexsurv.response(model.extract(m, "response"))
    mml <- mx <- vector(mode="list", length=length(dlist$pars))
    names(mml) <- names(mx) <- c(dlist$location, setdiff(dlist$pars, dlist$location))
    for (i in names(forms)){
        mml[[i]] <- model.matrix(forms[[i]], m)
        mx[[i]] <- length(unlist(mx)) + seq_len(ncol(mml[[i]][,-1,drop=FALSE]))
    }
    X <- flexsurv:::compress.model.matrices(mml)

    weights <- model.extract(m, "weights")
    if (is.null(weights)) weights <- m$"(weights)" <- rep(1, nrow(X))
    bhazard <- model.extract(m, "bhazard")
    if (is.null(bhazard)) bhazard <- rep(0, nrow(X))
    dat <- list(Y=Y, m=m, mml=mml)
    ncovs <- length(attr(m, "covnames.orig"))

    ncoveffs <- ncol(X)
    nbpars <- length(parnames) # number of baseline parameters
    npars <- nbpars + ncoveffs

    if (missing(inits) && is.null(dlist$inits))
        stop("\"inits\" not supplied, and no function to estimate them found in the custom distribution list")
    if (missing(inits) || any(is.na(inits))){
        yy <- ifelse(Y[,"status"]==3 & is.finite(Y[,"time2"]), (Y[,"time1"] + Y[,"time2"])/2, Y[,"time1"])
        wt <- yy*weights*length(yy)/sum(weights)
        dlist$inits <- flexsurv:::expand.inits.args(dlist$inits)
        inits.aux <- c(aux, list(forms=forms, data=if(missing(data)) NULL else data, weights=temp$weights,
                                 control=sr.control,
                                 counting=(attr(model.extract(m, "response"), "type")=="counting")
        ))
        auto.inits <- dlist$inits(t=wt,mf=m,mml=mml,aux=inits.aux)
        if (!missing(inits) && any(is.na(inits))) inits[is.na(inits)] <- flexsurv:::auto.inits[is.na(inits)]
        else inits <- auto.inits
    }
    if (!is.numeric(inits)) stop ("initial values must be a numeric vector")
    nin <- length(inits)
    if (nin < npars && ncoveffs > 0)
        inits <- c(inits, rep(0,length=npars-nin))
    else if (nin > npars){
        inits <- inits[1:npars]
        warning("Initial values are a vector length > ", npars, ": using only the first ", npars)
    }
    else if (nin < nbpars){
        stop("Initial values are a vector length ", nin, ", but distribution has ",nbpars, " parameters")
    }

    for (i in 1:nbpars)
        inits[i] <- dlist$transforms[[i]](inits[i])
    outofrange <- which(is.nan(inits) | is.infinite(inits))
    if (any(outofrange)){
        plural <- if(length(outofrange) > 1) "s" else ""
        stop("Initial value", plural, " for parameter", plural, " ",
             paste(outofrange,collapse=","), " out of range")
    }
    names(inits) <- c(parnames, colnames(X))

    if (!is.null(fixedpars) && !is.logical(fixedpars) &&
        (!is.numeric(fixedpars) || any(!(fixedpars %in% 1:npars)))){
        dots <- if(npars>2) "...," else ""
        stop("fixedpars must be TRUE/FALSE or a vector of indices in 1,",dots,npars)
    }

    if ((is.logical(fixedpars) && fixedpars==TRUE) ||
        (is.numeric(fixedpars) && identical(fixedpars, 1:npars))) {
        minusloglik <- flexsurv:::minusloglik.flexsurv(inits, Y=Y, X=X,
                                                       weights=weights, bhazard=bhazard,
                                                       dlist=dlist, inits=inits,
                                                       dfns=dfns, aux=aux, mx=mx)
        res.t <- matrix(inits, ncol=1)
        inits.nat <- inits
        for (i in 1:nbpars)
            inits.nat[i] <- dlist$inv.transforms[[i]](inits[i])
        res <- matrix(inits.nat, ncol=1)
        dimnames(res) <- dimnames(res.t) <- list(names(inits), "est")
        ret <- list(res=res, res.t=res.t, npars=0,
                    loglik=-as.vector(minusloglik), logliki=attr(minusloglik,"indiv"))
    }
    else {
        optpars <- inits[setdiff(1:npars, fixedpars)]
        optim.args <- list(...)
        if (is.null(optim.args$method)){
            optim.args$method <- "BFGS"
        }
        gr <- if (dfns$deriv) flexsurv:::Dminusloglik.flexsurv else NULL

        logLikFactory2(Y=Y, X=X,
                       weights=weights,
                       bhazard=bhazard,
                       inits=inits, dlist=dlist,
                       dfns=dfns,
                       aux=aux, mx=mx,
                       fixedpars=fixedpars)


    }

}




flexsurvspline.likelihood <- function(formula, data, weights, bhazard, subset,
                                      k=0, knots=NULL, bknots=NULL, scale="hazard", timescale="log", ...){
    ## Get response matrix from the formula.  Only need this to obtain
    ## default knots.  Largely copied from flexsurvreg - ideally
    ## should be in separate function, but can't make scoping work.

    call <- match.call()
    indx <- match(c("formula", "data", "weights", "bhazard", "subset", "na.action"), names(call), nomatch = 0)
    if (indx[1] == 0)
        stop("A \"formula\" argument is required")
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    ## remove the predictors
    f2 <- as.formula(gsub("(.*~).*", "\\1 1", Reduce(paste, deparse(formula))))
    environment(f2) <- environment(formula)
    temp[["formula"]] <- f2
    if (missing(data)) temp[["data"]] <- environment(formula)
    if (missing(data)) data <- environment(formula) # TESTME to pass to flexsurvreg
    m <- eval(temp, parent.frame())
    Y <- flexsurv:::check.flexsurv.response(model.extract(m, "response"))

    dtimes <- Y[,"stop"][Y[,"status"]==1]
    if (is.null(knots)) {
        is.wholenumber <-
            function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
        if (is.null(k)) stop("either \"knots\" or \"k\" must be specified")
        if (!is.numeric(k)) stop("k must be numeric")
        if (!is.wholenumber(k) || (k<0)) stop("number of knots \"k\" must be a non-negative integer")
        knots <- quantile(flexsurv:::tsfn(dtimes,timescale), seq(0, 1, length=k+2)[-c(1,k+2)])
    }
    else {
        if (!is.numeric(knots)) stop("\"knots\" must be a numeric vector")
        minlogtime <- min(tsfn(Y[,"stop"], timescale))
        if (any(knots <= minlogtime)) {
            badknots <- knots[knots < min(tsfn(Y[,"stop"],timescale))]
            plural <- if (length(badknots) > 1) "s" else ""
            stop(sprintf("knot%s %s less than or equal to minimum %stime", plural, paste(badknots,collapse=", "), (if (timescale=="log") "log " else "")))
        }
        maxlogtime <- max(tsfn(Y[,"stop"], timescale))
        if (any(knots >= maxlogtime)) {
            badknots <- knots[knots > maxlogtime]
            plural <- if (length(badknots) > 1) "s" else ""
            stop(sprintf("knot%s %s greater than or equal to maximum %stime", plural, paste(badknots,collapse=", "), (if (timescale=="log") "log " else "")))
        }
    }
    if (is.null(bknots)) {
        ## boundary knots chosen from min/max observed death times...
        if (length(dtimes)>0) {
            bt <- dtimes
        } else {
            ## unless all observations censored, where censoring times used instead
            ## "time" used with right censoring
            ## "time1", "time2" used with interval censoring
            bt <- c(Y[,"time1"], Y[,"time2"], Y[,"time"])
            bt <- bt[is.finite(bt)]
        }
        bknots <- c(min(flexsurv:::tsfn(bt,timescale)), max(flexsurv:::tsfn(bt,timescale)))
        if (bknots[1] == bknots[2])
            warning("minimum and maximum log death times are the same: knot and boundary knot locations should be supplied explicitly")
    } else
        if (!is.numeric(bknots) || (length(bknots) !=2) ) stop("boundary knots should be a numeric vector of length 2")
    knots <- c(bknots[1], knots, bknots[2])

    nk <- length(knots)
    custom.fss <- list(
        name = "survspline", # unused, d,p functions passed through
        pars = c(paste0("gamma",0:(nk-1))),
        location = c("gamma0"),
        transforms = rep(c(identity), nk), inv.transforms=rep(c(identity), nk),
        inits = flexsurv:::flexsurv.splineinits
    )
    aux <- list(knots=knots, scale=scale, timescale=timescale)
    dfn <- unroll.function(dsurvspline, gamma=0:(nk-1))
    pfn <- unroll.function(psurvspline, gamma=0:(nk-1))
    rfn <- unroll.function(rsurvspline, gamma=0:(nk-1))
    hfn <- unroll.function(hsurvspline, gamma=0:(nk-1))
    Hfn <- unroll.function(Hsurvspline, gamma=0:(nk-1))
    qfn <- unroll.function(qsurvspline, gamma=0:(nk-1))
    meanfn <- unroll.function(mean_survspline, gamma=0:(nk-1))
    rmstfn <- unroll.function(rmst_survspline, gamma=0:(nk-1))
    Ddfn <- if (scale=="normal") NULL else unroll.function(flexsurv:::DLdsurvspline, gamma=0:(nk-1))
    DSfn <- if (scale=="normal") NULL else unroll.function(flexsurv:::DLSsurvspline, gamma=0:(nk-1))
    args <- c(list(formula=formula, data=data, dist=custom.fss,
                   dfns=list(d=dfn,p=pfn,r=rfn,h=hfn,H=Hfn,rmst=rmstfn,mean=meanfn, q=qfn,
                             DLd=Ddfn,DLS=DSfn,deriv=!(scale=="normal")), aux=aux), list(...))


    ## Try an alternative initial value routine if the default one gives zero likelihood
    fpold <- args$fixedpars
    args$fixedpars <- TRUE
    #  if (is.infinite(do.call("flexsurvreg", args)$loglik)){
    #      args$dist$inits <- flexsurv:::flexsurv.splineinits.cox
    #  }
    args$fixedpars <- fpold
    args$weights <- temp$weights
    args$bhazard <- temp$bhazard
    args$subset <- temp$subset



    ret <- do.call("flexsurvreg.likelihood", args) # faff to make ... args work within functions

    ret
}



