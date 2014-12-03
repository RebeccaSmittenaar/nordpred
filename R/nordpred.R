
#' Fit power5 and poisson Age-Period-Cohort models for prediction of cancer incidence
#' 
#' \code{nordpred} uses the power5 and poisson Age-Period-Cohort (APC) models to 
#' calculate prediction of cancer incidence and mortality.
#' 
#' For details of the choice of prediction base, significance test for using 
#' recent slope, and for the power5 model, see Moller B., Fekjaer H. et al. (2002), 
#' see references.
#' 
#' @param cases A \code{data.frame} with number of cases
#' @param pyr A \code{data.frame} with observed and forecasted person years.
#' @param startestage Youngest age group to be included in the regression model. 
#' Predictions for age groups below this limit it based on average rates from 
#' the last 10 years.
#' @param startuseage Youngest age group which uses regression model as basis 
#' for predicted rates
#' @param noperiods A list of candidate number of periods in prediction base 
#' (e.g \code{4:6}).If the goodness of fit test is rejected based on the widest base 
#' (e.g. 6 periods), the first period is exclude etc. Use a fixed number to force 
#' a specific prediction base. If e.g. \code{noperiods = 5}, predictions is based on the 
#' last 5 five-year periods, irrespective of the result a goodness of fit evaluation
#' @param recent Project average trend or use the slope for the last 10 years? 
#' (If \code{recent = FALSE}, average trend for the whole observation period is used, if \code{recent = TRUE}, 
#' the slope from the last 10 years is used. If \code{NULL} (default) the choice is based 
#' on a significance test for departure from linear trend
#' @param cuttrend Cut trend in predictions? Default is 0 \%, 25 \%, 50 \%, 75 \%, 
#' 75 \% cut in drift (a vector of proportions of drift to cut in each projection period)
#' @param linkfunc Link function to use in the model. Default is special 
#' version used in the Nordpred project ('power5'), where the link is \eqn{g(x)=x^0.2}, 
#' while the alternative is the poisson function ('poisson'), where the link is 
#' \eqn{g(x)=log(x)}
#' @param x an object to test for class \code{nordpred}
#' 
#' @return 
#' 
#' \code{nordpred} returns an object of class \code{nordpred} 
#' (see \code{\link{nordpred.object}}). 
#' 
#' \code{is.nordpred} returns \code{TRUE} if input object is of class 
#' \code{nordpred}, \code{FALSE} otherwise.
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
#' # Fit model & predict new incidence:
#' res <- nordpred(indata, inpop, startestage = 5, startuseage = 6, 
#'          cuttrend = c(0, .25, .5, .75, .75))
#' res2 <- nordpred(indata, inpop, startestage = 5, startuseage = 6,
#'      cuttrend = c(0, .25, .5, .75, .75), linkfunc = 'poisson')
#'      
#' # Print / get results: 
#' print(res) 
#' nordpred.getpred(res) 
#' summary(res, printpred = FALSE)
#' 
#' # Get results with standardisation:
#' wstand <- c(0.12, 0.1, 0.09, 0.09, 0.08, 0.08, 0.06, 0.06, 0.06, 0.06,0.05,
#'             0.04, 0.04, 0.03, 0.02, 0.01, 0.005, 0.005) 
#' round(nordpred.getpred(res, incidence = TRUE, standpop = NULL), 2) 
#' round(nordpred.getpred(res, incidence = TRUE, standpop = wstand), 2)
#' 
#' # Plot results: 
#' plot(res, standpop = wstand)
#' 
#' # Plot results with power5 and poisson links: 
#' plot(res2, standpop = wstand) 
#' plot(res, new = FALSE, lty = c(1,2), standpop = wstand)
#' 
#' # Make estimates:
#' est <- nordpred.estimate(cases = indata, pyr = inpop, noperiod = 4, startestage = 5)
#' 
#' # Different cut trend scenarios, using average drift (recent = FALSE): 
#' plot(nordpred.prediction(est, startuseage = 6, cuttrend = c(0, 0, 0, 0, 0), recent = FALSE),
#'      standpop = wstand, new = TRUE) 
#' plot(nordpred.prediction(est, startuseage = 6, cuttrend = c(1, 1, 1, 1, 1), recent = FALSE),
#'      standpop = wstand, new = FALSE, lty = c(1, 2)) 
#' plot(nordpred.prediction(est, startuseage = 6, cuttrend = c(0, .25, .5, .75, .75),
#'      recent = FALSSE), standpop = wstand, new = FALSE, lty = c(1, 4))
#'      
#' @export
#' @family nordpred
#' @name nordpred


nordpred <- function(cases, pyr, startestage, startuseage, noperiods = NULL, recent = NULL, 
    cuttrend = c(0, 0.25, 0.5, 0.75, 0.75), linkfunc = "power5") {
    
    Rplatform <- exists("R.Version")
    
    # Seting default and checking data:
    percases <- dim(cases)[2]
    # Number of periods in prediction base (observed periods)
    
    if (percases < 3) {
        stop("Too few periods in \"cases\"")
    }
    
    if (is.null(noperiods)) {
        noperiods <- c(min(percases, 4):min(percases, 6))
        # List possible candidates for number of periods to base predictions on. Default
        # is 4:6 if available
    }
    
    # Choose number of periods by cutting stepwise execution of the highest candidate
    # number (i.e. cutting the most ancient periods)
    noperiods <- sort(noperiods)
    while (length(noperiods) > 1) {
        maxnoperiod <- max(noperiods)
        glm <- nordpred.estimate(cases, pyr, maxnoperiod, startestage)$glm
        if (Rplatform) {
            pvalue <- 1 - pchisq(glm$deviance, glm$df.residual)
        } else {
            pvalue <- 1 - pchisq(glm$deviance, glm$df)
        }
        if (pvalue < 0.01) {
            noperiods <- noperiods[1:(length(noperiods) - 1)]
        } else {
            noperiods <- maxnoperiod
        }
    }
    noperiod <- noperiods
    
    # Set status for recent (whether to use recent trend or average trend)
    if (is.null(recent)) {
        recent <- nordpred.estimate(cases, pyr, noperiod, startestage)$suggestionrecent
    }
    
    # Perform estimation and prediction:
    est <- nordpred.estimate(cases = cases, pyr = pyr, noperiod = noperiod, startestage = startestage, 
        linkfunc = linkfunc)
    pred <- nordpred.prediction(nordpred.estimate.object = est, startuseage = startuseage, 
        recent = recent, cuttrend = cuttrend)
    return(pred)
} 


#' @rdname nordpred
#' @export
is.nordpred <- function(x){
    inherits(x, "nordpred")
}
