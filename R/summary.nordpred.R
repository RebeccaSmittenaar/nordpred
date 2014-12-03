

#' Makes a summary of a nordpred object
#' 
#' \code{summary.nordpred} uses a \code{nordpred} object (see \code{\link{nordpred.object}})
#' to summarize 
#' the information
#' 
#' 
#' @param nordpred.object An object of class \code{nordpred} (see \code{\link{nordpred.object}}) 
#' @param printpred Indicates whether to print the observed and predicted number of cases
#' @param printcall Indicates whether to print the function call
#' @param digits Specifies the number of digits in the tabulation
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
#' res <- nordpred(indata, inpop, startestage = 5, startuseage = 6, cuttrend = c(0,.25,.5,.75,.75))
#' res2 <- nordpred(indata, inpop, startestage = 5, startuseage = 6,
#'      cuttrend = c(0, .25, .5, .75, .75), linkfunc = 'poisson')
#'      
#' # Print / get results: 
#' print(res) 
#' nordpred.getpred(res) 
#' summary(res, printpred = FALSE)
#' 
#' @export
#' @family nordpred


summary.nordpred <- function(nordpred.object, printpred = TRUE, printcall = FALSE, digits = 1) {
    
    if (class(nordpred.object) != "nordpred") {
        stop("Variable \"nordpred.object\" must be of type \"nordpred\"")
    }
    
    # Setting internal variables:
    obsto <- names(nordpred.object$predictions)[dim(nordpred.object$predictions)[2] - 
        nordpred.object$nopred]
    
    if (!is.null(nordpred.object$pvaluerecent)) {
        precent <- round(nordpred.object$pvaluerecent, 4)
    } else {
        precent <- NA
    }
    
    if (!is.null(nordpred.object$gofpvalue)) {
        gofpvalue <- round(nordpred.object$gofpvalue, 4)
    } else {
        gofpvalue <- NA
    }
    
    # Print information about object:
    if (printpred) {
        cat("Observed and predicted values:")
        cat("(observations up to", obsto, ")\n")
        print(round(as.matrix(nordpred.object$predictions), digits = digits))
        cat("\n")
    }
    cat("\nPrediction done with:\n")
    
    moptions <- matrix(NA, 8, 2)
    moptions[, 1] <- c("Number of periods predicted (nopred):", "Trend used in predictions (cuttrend):", 
        "Number of periods used in estimate (noperiod):", "P-value for goodness of fit:", 
        "Used recent (recent):", "P-value for recent:", "First age group used (startuseage):", 
        "First age group estimated (startestage):")
    moptions[, 2] <- c(nordpred.object$nopred, paste(nordpred.object$cuttrend, collapse = " , "), 
        nordpred.object$noperiod, gofpvalue, nordpred.object$recent, precent, nordpred.object$startuseage, 
        nordpred.object$startestage)
    maxl <- max(nchar(moptions[, 1]))
    
    for (i in 1:dim(moptions)[1]) {
        spaces <- rep(" ", maxl - nchar(moptions[i, 1]) + 2)
        cat(moptions[i, 1], spaces, moptions[i, 2], "\n", sep = "")
    }
    
    if (printcall) {
        cat("\n  Call: ")
        dput(attr(nordpred.object, "Call"))
    }
    invisible(nordpred.object)
}
 
