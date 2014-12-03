
#' Calculates predictions based on a nordpred.estimate object
#'  
#'  \code{nordpred.prediction} uses a \code{\link{nordpred.estimate.object}} 
#'  to calculate prediction of cancer incidence and mortality
#' 
#'@param nordpred.estimate.object A \code{\link{glm}}-object based on 
#' \code{\link{nordpred.estimate}} (see \code{\link{nordpred.estimate.object}}).
#' @param startuseage Youngest age group which uses regression model as basis 
#' for predicted rates 
#' @param recent Project average trend or use the slope for the last 10 years? 
#' (If \code{recent = FALSE}, average trend for the whole observation period is
#' used, if \code{recent = TRUE}, the slope from the last 10 years is used)
#' @param cuttrend Cut trend in predictions? (a vector of proportionsof drift
#' to be cut in each projection period)
#' 
#' @return an object of class \code{nordpred} (see \code{\link{nordpred.object}}).
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
#' # Use estimat object to make predictions:
#' res <- nordpred.prediction(est,startuseage = 6, cuttrend = c(0, .25, .5, .75, .75), recent = TRUE)
#' res2 <- nordpred.prediction(est2, startuseage = 6, cuttrend = c(0, .25, .5, .75, .75), recent = TRUE)
#' 
#' # Get results: 
#' print(res) 
#' nordpred.getpred(res) 
#' summary(res, printpred = FALSE)
#' 
#' @export
#' @family nordpred


nordpred.prediction <- function(nordpred.estimate.object, startuseage, recent, 
                                cuttrend = c(0, 0.25, 0.5, 0.75, 0.75)) {
    
    if (class(nordpred.estimate.object) != "nordpred.estimate") {
        stop("Variable \"nordpred.estimate.object\" must be of type \"nordpred.estimate\"")
    }
    
    if (nordpred.estimate.object$startestage > startuseage) {
        stop("\"startuseage\" is set to high compared to \"startestage\" in \"nordpred.estimate.object\"")
    }
    
    if (length(cuttrend) < (dim(nordpred.estimate.object$pyr)[2] - dim(nordpred.estimate.object$cases)[2])) {
        err <- paste("\"cuttrend\" must always be at least the same length as")
        err <- paste(err, "the number of periods with population forecasts")
        stop(err)
    }
    
    # Setting local variables:
    cases <- nordpred.estimate.object$cases
    pyr <- nordpred.estimate.object$pyr
    noperiod <- nordpred.estimate.object$noperiod
    nototper <- dim(pyr)[2]
    noobsper <- dim(cases)[2]
    nonewpred <- nototper - noobsper
    
    if (length(cuttrend) > nonewpred) {
        cuttrend <- cuttrend[1:nonewpred]
    }
    
    if (is.data.frame(pyr)) {
        years <- names(pyr)
    } else {
        if (is.null(dimnames(pyr))) {
            years <- paste("Periode", 1:nototper)
        } else {
            years <- dimnames(pyr)[[2]]
        }
    }
    
    # We make data object:
    datatable <- matrix(NA, 18, nototper)
    datatable[, 1:(nototper - nonewpred)] <- as.matrix(cases)
    datatable <- data.frame(datatable)
    row.names(datatable) <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", 
        "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", 
        "80-84", "85+")
    names(datatable) <- years
    # We fill in observed cases:
    
    # Calculate predictions in number of cases:
    for (age in 1:(startuseage - 1)) {
        # For agegroups with litle data, we use mean incidence for last two periods:
        obsinc <- cases[age, (noobsper - 1):noobsper]/pyr[age, (noobsper - 1):noobsper]
        if (sum(is.na(obsinc))) {
            obsinc[is.na(obsinc)] <- 0
        }
        datatable[age, (noobsper + 1):nototper] <- ((obsinc[, 1] + obsinc[, 2])/2) * 
            pyr[age, (noobsper + 1):nototper]
    }
    for (age in startuseage:18) {
        startestage <- nordpred.estimate.object$startestage
        coefficients <- nordpred.estimate.object$glm$coefficients
        
        coh <- (18 - startestage) - (age - startestage) + (noperiod + 1:nonewpred)
        # No. agegoups - age + period
        noages <- 18 - startestage + 1
        driftmp <- cumsum(1 - cuttrend)
        cohfind <- noages + (noperiod - 1) + 1 + (coh - 1)
        maxcoh <- 18 - startuseage + noperiod
        agepar <- as.numeric(coefficients[age - startestage + 1])
        driftfind <- pmatch("Period", attributes(coefficients)$names)
        driftpar <- as.numeric(coefficients[driftfind])
        
        cohpar <- rep(NA, length(coh))
        for (i in 1:length(coh)) {
            if (coh[i] < maxcoh) {
                cohpar[i] <- as.numeric(coefficients[cohfind[i]])
            } else {
                # For young cohorts not estimated:
                cohpar[i] <- as.numeric(coefficients[length(coefficients) - (startuseage - 
                  startestage)])
                cohpar[i][is.na(cohpar[i])] <- 0
            }
        }
        
        if (recent) {
            lpfind <- driftfind + noperiod - 2
            # -2 because first and last period effect is set as contrast (p.first=p.last)
            lppar <- as.numeric(coefficients[lpfind])
            # lppar is the slope corresponding to a linear trend in the last 10 years (two
            # periods)
            driftrecent <- driftpar - lppar
        }
        
        if (nordpred.estimate.object$linkfunc == "power5") {
            if (recent) {
                rate <- (agepar + driftpar * noobsper + driftrecent * driftmp + cohpar)^5
            } else {
                rate <- (agepar + driftpar * (noobsper + driftmp) + cohpar)^5
            }
        } else {
            # Possion link:
            if (recent) {
                rate <- exp(agepar + driftpar * noobsper + driftrecent * driftmp + 
                  cohpar)
            } else {
                rate <- exp(agepar + driftpar * (noobsper + driftmp) + cohpar)
            }
        }
        datatable[age, (noobsper + 1):nototper] <- rate * pyr[age, (noobsper + 1):nototper]
    }
    
    # Structure and return results:
    res <- list(predictions = datatable, pyr = pyr, nopred = nonewpred, noperiod = nordpred.estimate.object$noperiod, 
        gofpvalue = nordpred.estimate.object$gofpvalue, recent = recent, pvaluerecent = nordpred.estimate.object$pvaluerecent, 
        cuttrend = cuttrend, startuseage = startuseage, startestage = nordpred.estimate.object$startestage, 
        glm = nordpred.estimate.object$glm)
    class(res) <- "nordpred"
    attr(res, "Call") <- sys.call()
    return(res)
} 
