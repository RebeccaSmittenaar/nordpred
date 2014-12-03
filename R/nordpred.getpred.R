
#' Gets the observed and predicted incidence rates on matrix form
#'  
#' \code{nordpred.getpred} uses a nordpred object to extract the observed 
#' and predicted incidence rates
#' 
#' @param nordpred.object An object based on \code{\link{nordpred}} or
#'  \code{\link{nordpred.prediction}}.
#' @param incidence Indicates whether to give incidence or number of cases
#' @param standpop A vector of weights for age standardisation. 
#' Default is no standardisation (crude rates), but using a standardisation 
#' (for the suitable no of age groups) is recommended
#' @param excludeobs Exclude number for observed periods and only give numbers 
#' for predicted periods
#' @param byage Report numbers by age groups. If false, crude or age 
#' standardised rates are given
#' @param agegroups Which agegroups to include. E.g. \code{c(5:18)} 
#' includes age groups five to eighteen
#' 
#' @return an object of class \code{\link{nordpred}}.
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
#' # Reading package: 
#' source("nordpred.S")
#' 
#' # Reading data (Colon cancer for Norwegian males)
#' indata <- read.table("data//colon-men-Norway.txt",header =T,sep=",",row.names=1)
#' inpop1 <- read.table("data//men-Norway.txt",header =T,sep=",",row.names=1)
#' inpop2 <- read.table("data//men-Norway-pred.txt",header =T,sep=",",row.names=1)
#' 
#' # Include possible population predictions 
#' inpop <- cbind(inpop1,inpop2)
#' 
#' # Fit model & predict new incidence:
#' res <- nordpred(indata,inpop,startestage=5,startuseage=6,cuttrend=c(0,.25,.5,. 75,.75))
#' res2 <- nordpred(indata,inpop,startestage=5,startuseage=6,cuttrend=c(0,.25,.5,. 75,.75),linkfunc="poisson")
#' 
#' # Print / get results: 
#' print(res) 
#' nordpred.getpred(res) 
#' summary(res,printpred=F)
#' 
#' # Get results with standardisation:'
#' wstand <- c(0.12, 0.1, 0.09, 0.09, 0.08, 0.08, 0.06, 0.06, 0.06, 0.06,0.05,
#'             0.04, 0.04, 0.03, 0.02, 0.01, 0.005, 0.005) 
#' round(nordpred.getpred(res,incidence=T,standpop=NULL),2) 
#' round(nordpred.getpred(res,incidence=T,standpop=wstand),2)
#' 
#' @export
#' @family nordpred


nordpred.getpred <- function(nordpred.object,incidence=T,standpop=NULL,excludeobs=F,byage,agegroups="all") {
       
    # Seting defaults:
    if (missing(byage)) {
        byage <- ifelse(is.null(standpop),T,F)
    }
    
    # Checking imput:
    if (class(nordpred.object)!="nordpred") {
        stop("Variable \"nordpred.object\" must be of type \"nordpred\"")	
    } 
    
    if ((!is.null(standpop)) && (!incidence)) {
        stop("\"standpop\" should only be used with incidence predictions (incidence=T)") 	
    }
    
    if (!is.null(standpop)) {
        if (round(sum(standpop),5)!=1) {
            stop("\"standpop\" must be of sum 1") 	
        }
        if ((length(standpop)!=length(agegroups)) && (agegroups[1]!="all")) {
            stop("\"standpop\" must be the same length as \"agegroups\"") 	
        }
        if (byage) {
            stop("\"standpop\" is only valid for \"byage=T\"") 
        }
    }
    
    # Seting local data:
    datatable <- nordpred.object$predictions
    pyr       <- data.frame(nordpred.object$pyr)
    
    # Secting agegroups:
    if (agegroups[1]!="all") {
        datatable <- datatable[agegroups,]
        pyr       <- pyr[agegroups,]
    }
    
    # If needed; Standardize data and Collapse agegroups
    if (!is.null(standpop)) {
        datainc <- (datatable/pyr)*100000
        if (sum(is.na(datainc))>0) {
            datainc[is.na(datainc)] <- 0
        } 	
        res <- apply(datainc*standpop,2,sum) 
    } else {
        if (!byage) {
            datatable <- apply(datatable,2,sum)
            pyr       <- apply(pyr,2,sum) 
        }
        if (incidence) {
            res <- (datatable/pyr)*100000
            if (sum(is.na(res))>0) {
                res[is.na(res)] <- 0
            } 	
        } else {
            res <- datatable
        } 
    }
    
    # Select data:
    if (excludeobs) {
        if (is.matrix(res)) {
            predstart <- dim(res)[2]-nordpred.object$nopred+1
            res <- res[,predstart:(predstart+nordpred.object$nopred-1)]
        } else {
            predstart <- length(res)-nordpred.object$nopred+1
            res <- res[predstart:(predstart+nordpred.object$nopred-1)]
        }
    } 
    
    # Return data: 
    return(res)
}
