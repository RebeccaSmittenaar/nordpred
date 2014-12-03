
#' Prints a nordpred object
#' 
#' \code{print.nordpred} prints the observed and predicted number of cases in 
#' a nordpred object
#' 
#' @param nordpred.object An object of class \code{nordpred} (see \code{\link{nordpred.object}})
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
#' res <- nordpred(indata, inpop, startestage = 5, startuseage = 6,
#'      cuttrend = c(0, .25, .5, .75, .75))
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


print.nordpred <- function(nordpred.object, digits = 1) {
    
    if (class(nordpred.object) != "nordpred") {
        stop("Variable \"nordpred.object\" must be of type \"nordpred\"")
    }
    
    # Setting internal variables:
    obsto <- names(nordpred.object$predictions)[dim(nordpred.object$predictions)[2] - 
        nordpred.object$nopred]
    # Print information about object:
    cat("Observed and predicted values:")
    cat("(observations up to", obsto, ")\n")
    print(round(as.matrix(nordpred.object$predictions), digits = digits))
    cat("\n  Call: ")
    dput(attr(nordpred.object, "Call"))
    invisible(nordpred.object)
    
} 
