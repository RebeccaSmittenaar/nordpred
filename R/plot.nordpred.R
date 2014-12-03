
#' plots the predicted rates from a nordpred object
#' 
#' \code{plot.nordpred} uses nordpred object to plot observed and predicted rates
#' 
#' This function is a method for the generic function plot for class nordpred. 
#' It can be invoked by calling \code{\link{plot}} for an object x of the appropriate class, 
#' or directly by calling \code{plot.nordpred} regardless of the class of the object. 
#' For more available options, see \code{\link{plot}}. 
#' For details of the choice of prediction base, significance test for using 
#' recent slope, and for the power5 model, see Moller B., Fekjaer H. et al. (2002), 
#' see references.
#' 
#' @param nordpred.object An object of class \code{nordpred} (see \code{\link{nordpred.object}}
#' @param incidence Indicates whether to plot incidence or number of cases
#' @param standpop A vector of weights for age standardisation. 
#' Default is no standardisation (crude rates), but using a standardisation 
#' (for the suitable no of age groups) is recommended
#' @param agegroups Which agegroups to include
#' @param noperiods A list of candidate number of periods in prediction base 
#' (e.g \code{4:6}).If the goodness of fit test is rejected based on the widest base 
#' (e.g.6 periods), the first period is exclude etc. Use a fixed number to force 
#' a specific prediction base. If e.g. \code{noperiods = 5}, predictions is based on the 
#' last 5 five-year periods, irrespective of the result a goodness of fit evaluation
#' @param recent Project average trend or use the slope for the last 10 years? 
#' (If \code{recent = FALSE}, average trend for the whole observation period is used, if \code{recent = TRUE}, 
#' the slope from the last 10 years is used. If \code{NULL} (default) the choice is based 
#' on a significance test for departure from linear trend
#' @param cuttrend Cut trend in predictions? Default is 0 \%, 25 \%, 50 \%, 75 \%, 
#' 75 \% cut in drift (a vector of proportions of drift to cut in each projection period)
#' @param linkfunc Link function to use in the model. Default is special 
#' version used in the Nordpred project ('power5'), where the link is \eqn{g(x) = x^0.2}, 
#' while the alternative is the poisson function ('poisson'), where the link is 
#' \eqn{g(x) = log(x)}
#' 
#' @return object of class \code{nordpred} (see \code{\link{nordpred.object}}).
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
#' # Get results with stanardiziotion:
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
#' plot(res, new = FALSE, lty = c(1, 2), standpop = wstand)
#' 
#' # Different cut trend scenarios, using average drift (recent = FALSE): 
#' plot(nordpred.prediction(est, startuseage = 6, cuttrend = c(0,0,0,0,0), 
#'      recent = FALSE), standpop = wstand, new = TRUE) 
#' plot(nordpred.prediction(est, startuseage = 6, cuttrend = c(1,1,1,1,1),
#'      recent = FALSE), standpop = wstand, new = FALSE, lty = c(1,2)) 
#' plot(nordpred.prediction(est,startuseage=6,cuttrend=c(0,.25,.5,.75,.75),
#'      recent = FALSE), standpop = wstand, new = FALSE, lty = c(1,4))      
#'      
#' @export
#' @family nordpred


plot.nordpred <- function(nordpred.object, incidence = TRUE, standpop = NULL, agegroups = "all", 
    startplot = 1, xlab = "", ylab = "", main = "", labels = NULL, ylim = NULL, lty = c(1, 
        3), col = c(1, 1), new = TRUE, ...) {
    
    if (class(nordpred.object) != "nordpred") {
        stop("Variable \"nordpred.object\" must be of type \"nordpred\"")
    }
    
    # Seting internal variables:
    nopred <- nordpred.object$nopred
    if (is.null(labels)) {
        labels <- dimnames(nordpred.object$predictions)[[2]]
        labels <- labels[startplot:length(labels)]
    }
    
    # Reding & formating data:
    indata <- nordpred.getpred(nordpred.object, incidence = incidence, standpop = standpop, 
        agegroups = agegroups, byage = F)
    indata <- indata[startplot:length(indata)]
    
    # Create plots:
    maxx <- length(indata)
    if (new) {
        maxy <- max(indata)
        if (is.null(ylim)) {
            ylim <- c(0, maxy)
        }
        plot(c(1, maxx), ylim, type = "n", ylab = ylab, xlab = xlab, axes = F, ...)
        axis(2)
        axis(1, at = 1:maxx, labels = labels)
        box()
        title(main)
    }
    lines(1:(maxx - nopred), indata[1:(maxx - nopred)], lty = lty[1], col = col[1], 
        ...)
    lines((maxx - nopred):maxx, indata[(maxx - nopred):maxx], lty = lty[2], col = col[2], 
        ...)
    
    # Returning object as invisible
    invisible(nordpred.object)
} 
