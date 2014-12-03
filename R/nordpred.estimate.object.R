
#' Nordpred.estimate-object with fit of power5 and poisson Age-Period-Cohort 
#' models for prediction of cancer incidence
#'  
#' \code{nordpred} uses the power5 and poisson Age-Period-Cohort (APC) models to 
#' calculate prediction of cancer incidence and mortality 
#' This class of objects is returned by the nordpred.estimate class of functions 
#' to represent a fit of power5 and poisson Age-Period-Cohort models for 
#' prediction of cancer incidence. 
#' Objects of this class have methods \code{\link{print.nordpred}}, 
#' \code{\link{summary}} and \code{\link{plot}}.
#' 
#' @section Components:
#' \describe{
#' \item{glm}{Fitted glm-object}
#' \item{cases}{A \code{data.frame} with number of cases}
#' \item{pyr}{A \code{data.frame} with observed and forecasted person years}
#' \item{noperiod}{Number of periods used in estimate}
#' \item{gofpvalue}{P-value for goodness of fit}
#' \item{startestage}{Youngest age group which have been included in the 
#' regression model. Predictions for age groups below this limit it based on 
#' average rates from the last 10 years.}
#' \item{suggestionrecent}{Indicator recommendation build on pvaluerecent for 
#' projecting of average trend or use the slope for the last 10 years? 
#' If recent=F, recommendation is to use average trend for the whole observation 
#' period, and if recent=T recommendation is to use the slope from the last 10 years}
#' \item{pvaluerecent}{P-value for use of recent trend based on a significance 
#' test for departure from linear trend}
#' \item{linkfunc}{Link function used in the model. 
#' Default is special version used in the Nordpred project ("power5"), 
#' where the link is \eqn{g(x)=x^0.2}, while the alternative is the poisson function 
#' ("poisson"), where the link is \eqn{g(x)=log(x)}}
#' }
#' 
#' The object will also contain the following (see \code{\link{lm}}): 
#' \code{formula}, \code{terms}, \code{assign} and \code{call}.
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
#' @family nordpred
#' @name nordpred.estimate.object
NULL

