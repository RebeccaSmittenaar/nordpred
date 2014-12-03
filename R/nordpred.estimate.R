
#' Fit power5 and poisson Age-Period-Cohort models for prediction of cancer incidence
#'  
#' \code{nordpred.estimate} estimates parameters in the power5 or poisson 
#' Age-Period-Cohort (APC) model
#' 
#' For details of the choice of prediction base, significance test for using 
#' recent slope, and for the power5 model, see Moller B., Fekjaer H. et al. (2002), 
#' see references.
#' 
#' @param cases A \code{data.frame} with number of cases
#' @param pyr A \code{data.frame} with observed and forecasted person years. 
#' @param startestage Youngest age group to include in the regression model
#' @param noperiod The number of periods to be used in prediction base.
#' @param lincfunc Link function to use in the model. 
#' Default is special version used in the Nordpred project ('power5'), 
#' where the link is \eqn{g(x) = x^0.2}, while the alternative is the poisson function 
#' ('poisson'), where the link is \eqn{g(x) = log(x)}
#' @param x an object to test for class \code{nordpred}

#' 
#' @return 
#' \code{nordpred.estimate} returns an object of class \code{nordpred.estimate} 
#' (see \code{\link{nordpred.estimate.object}}). 
#' 
#' \code{is.nordpred.estimate} returns \code{TRUE} if input object is of class 
#' \code{nordpred.estimate}, \code{FALSE} otherwise.
#' 
#' @references 
#' \itemize{
#' \item A website for nordpred is available at: 
#' \url{http://www.kreftregisteret.no/software/nordpred/}
#' \item Background for the methods can be found in: Moller B., Fekjaer H., Hakulinen T., 
#' Sigvaldason H, Storm H. H., Talback M. and Haldorsen T 'Prediction of cancer 
#' incidence in the Nordic countries: Empirical comparison of different approaches' 
#' Statistics in Medicine 2003; 22:2751-2766
#' \item An application of the function, using all the default settings, can be 
#' found in: Moller B, Fekjaer H, Hakulinen T, Tryggvadottir L, Storm HH, Talback M, 
#' Haldorsen T. Prediction of cancer incidence in the Nordic countries up to the 
#' year 2020. Eur J Cancer Prev Suppl 2002; 11: S1-S96
#' }
#' 
#' @author Harald Fekjaer and Bjorn Moller (Cancer Registry of Norway)
#' 
#' @section Note for S-plus:
#' Powerlink is made via a special modification in S-PLUS. This works fine 
#' for the point estimates, but the variance estimates found via the glm-objects 
#' are wrong. For variance estimates, we would rather recommend using R.
#' 
#' @examples
#' 
#' # data (Colon cancer for Norwegian males)
#' indata 
#' 
#' # Create dataset with observed and predicted population
#' inpop <- cbind(inpop1, inpop2)
#' 
#' # Fit model using powerlink (default):
#' est <- nordpred.estimate(cases = indata, pyr = inpop, noperiod = 4, startestage = 5)
#' 
#' # Fit model using poisson link:
#' est2 <- nordpred.estimate(indata, inpop, 4, 5, linkfunc = 'poisson')
#' 
#' # Print results: 
#' print(est) 
#' print(est$glm)
#' 
#' #' # Use estimat object to make predictions:
#' res <- nordpred.prediction(est, startuseage = 6, cuttrend = c(0, .25, .5, .75, .75), recent = TRUE)
#' 
#' @export
#' @family nordpred
#' @name nordpred.estimate


nordpred.estimate <- function(cases, pyr, noperiod, startestage, linkfunc = "power5", 
    find_best_model = FALSE) {
    
    Rplatform <- exists("R.Version")
    
    if (dim(cases)[1] != 18 || dim(pyr)[1] != 18) {
        stop("\"cases\" and \"pyr\" must have data for 18 age groups")
    }
    
    if (dim(cases)[2] > dim(pyr)[2]) {
        stop("\"pyr\" must include information about all periods in \"cases\"")
    }
    
    if (dim(pyr)[2] == dim(cases)[2]) {
        stop("\"pyr\" must include information on future rates (E.g. must be greater than \"cases\")")
    }
    
    if ((dim(pyr)[2] - dim(cases)[2]) > 5) {
        stop("Package can not predict more than 5 periods (given by sizes of \"pyr\" and \"cases\")")
    }
    
    if ((dim(cases)[2] - noperiod) < 0) {
        stop("More periods specified in \"noperiod\" than columns in \"cases\"")
    }
    
    if (noperiod < 3) {
        stop("\"noperiod\" must be 3 or larger to get enough for estimating")
    }
    
    # Setting internal variables:
    dnoperiods <- dim(cases)[2]
    dnoagegr <- 18
    
    # Transform dataformat:
    ageno <- rep(1:dnoagegr, dnoperiods)
    periodno <- sort(rep(1:dnoperiods, dnoagegr))
    cohort <- max(ageno) - ageno + periodno
    y <- c(as.matrix(pyr[, 1:dnoperiods]))
    apcdata <- data.frame(Cases = c(as.matrix(cases)), Age = ageno, Cohort = cohort, 
        Period = periodno, y = y)
    
    # Selecting data for regression:
    apcdata <- apcdata[apcdata$Age >= startestage, ]
    apcdata <- apcdata[apcdata$Period > (dnoperiods - noperiod), ]
    
    
    # Creation of power link:
    if (Rplatform) {
        # Sett variable for use in estimation
        y <- apcdata$y
        # Make link function:
        power5link <- poisson()
        power5link$link <- "0.2 root link Poisson family"
        power5link$linkfun <- function(mu) {
            (mu/y)^0.2
        }
        power5link$linkinv <- function(eta) {
            pmax(.Machine$double.eps, y * eta^5)
        }
        power5link$mu.eta <- function(eta) {
            pmax(.Machine$double.eps, 5 * y * eta^4)
        }
    } else {
        # Sett variable for use in estimation
        y <<- apcdata$y
        # Must be sett on top level Make link function:
        power5link <- poisson()
        power5link$link <- function(mu) {
            (mu/y)^0.2
        }
        power5link$inverse <- function(eta) {
            y * eta^(1/0.2)
        }
        power5link$deriv <- function(mu) {
            0.2 * (1/y) * (mu/y)^(0.2 - 1)
        }
    }
    
    # Setting contrast:
    options(contrasts = c("contr.treatment", "contr.poly"))
    
    # Estimation:
    if (linkfunc == "power5") {
        res.glm <- glm(Cases ~ as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) - 
            1, family = power5link, data = apcdata)
        if (find_best_model) {
            res.glm <- suppressWarnings(step(res.glm, trace = FALSE))
        }
    } else if (linkfunc == "poisson") {
        res.glm <- glm(Cases ~ as.factor(Age) + Period + as.factor(Period) + as.factor(Cohort) + 
            offset(log(y)) - 1, family = poisson(), data = apcdata)
    } else {
        stop("Unknown \"linkfunc\"")
    }
    
    if (Rplatform) {
        pvalue <- 1 - pchisq(res.glm$deviance, res.glm$df.residual)
    } else {
        pvalue <- 1 - pchisq(res.glm$deviance, res.glm$df)
    }
    
    # Basis for setting suggestion for 'recent' (whether to use recent trend or
    # average trend)
    mod1 <- glm(Cases ~ as.factor(Age) + Period + as.factor(Cohort) + offset(log(y)) - 
        1, family = poisson, data = apcdata)
    mod2 <- glm(Cases ~ as.factor(Age) + Period + I(Period^2) + as.factor(Cohort) + 
        offset(log(y)) - 1, family = poisson, data = apcdata)
    if (Rplatform) {
        pdiff <- anova(mod1, mod2, test = "Chisq")$"P(>|Chi|)"[2]
        # Correction added 2012-09-19 for compabilty with newer R versions:
        if (is.null(pdiff)) 
            pdiff <- anova(mod1, mod2, test = "Chisq")$"Pr(>Chi)"[2]
    } else {
        pdiff <- anova(mod1, mod2, test = "Chisq")[2, 7]
    }
    # use last two periods (5years*2) to predict if period square is significant,
    # else use whole drift
    if (pdiff < 0.05) {
        suggestionrecent <- T
    } else {
        suggestionrecent <- F
    }
    
    # Set class and return results
    res <- list(glm = res.glm, cases = cases, pyr = pyr, noperiod = noperiod, gofpvalue = pvalue, 
        startestage = startestage, suggestionrecent = suggestionrecent, pvaluerecent = pdiff, 
        linkfunc = linkfunc)
    class(res) <- "nordpred.estimate"
    attr(res, "Call") <- sys.call()
    return(res)
} 



#' @rdname nordpred.estimate
#' @export
is.nordpred.estimate <- function(x){
    inherits(x, "nordpred.estimate")
}
