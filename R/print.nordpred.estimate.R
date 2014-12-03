
#' Prints a nordpred estimate object
#' 
#' \code{print.nordpred.estimate} prints the estimation information from a 
#' \code{\link{nordpred.estimate.object}}.
#' 
#' @param nordpred.estimate.object An object of class \code{\link{nordpred.estimate.object}} 
#' 
#' @return object of class \code{\link{nordpred}}.
#' 
#' @references 
#' \itemize{
#' \item A website for nordpred is available at: 
#' \url{http://www.kreftregisteret.no/software/nordpred/}
#' \item Background for the methods can be found in: Moller B., Fekjaer H., Hakulinen T., 
#' Sigvaldason H, Storm H. H., Talback M. and Haldorsen T "Prediction of cancer 
#' incidence in the Nordic countries: Empirical comparison of different approaches" 
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
#' # Reading data (Colon cancer for Norwegian males)
#' indata <- read.table("data//colon-men-Norway.txt",header =T,sep=",",row.names=1)
#' inpop1 <- read.table("data//men-Norway.txt",header =T,sep=",",row.names=1)
#' inpop2 <- read.table("data//men-Norway-pred.txt",header =T,sep=",",row.names=1)
#' 
#' # Include possible population predictions 
#' inpop <- cbind(inpop1,inpop2)
#' 
#' # Fit model using powerlink (default):
#' est <- nordpred.estimate(cases=indata,pyr=inpop,noperiod=4,startestage=5)
#' 
#' # Fit model using poisson link:
#' est2 <- nordpred.estimate(indata,inpop,4,5,linkfunc="poisson")
#' 
#' # Print results: 
#' print(est) 
#' print(est$glm)
#' 
#' # Use estimat object to make predictions:
#' res <- nordpred.prediction(est,startuseage=6,cuttrend=c(0,.25,.5,.75,.75),rece nt=T)
#'      
#' @export
#' @family nordpred

print.nordpred.estimate <- function(nordpred.estimate.object) {
      
    if (class(nordpred.estimate.object)!="nordpred.estimate") {
        stop("Variable \"nordpred.estimate.object\" must be of type \"nordpred.estimate\"")	
    } 
    
    # Print information about object: 
    cat("Fit of cancer prediction model (as of Nordpred project)\n")
    cat("Fitted for",nordpred.estimate.object$noperiod,"periods")
    cat(", with estimation from age group number",nordpred.estimate.object$startestage,"\n")
    cat("  Call: ")
    dput(attr(nordpred.estimate.object,"Call"))
    
    invisible(nordpred.estimate.object)
}
