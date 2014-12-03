#' Incidence data for colon cancer among males in Norway
#'
#' A dataset containing observed number of colon cancer cases among Norweigen men.
#'
#' @format A data frame with 18 rows and 8 variables:
#' \describe{
#'   \item{rows}{Each row identifies a five year age interval, from 0-4 to 85+}
#'   \item{columns}{Each column identifies a five year period, from 1958-1962 to 1993-1997}
#' }
#' @examples
#' indata
#' @family nordpred_example_data
"indata"



#' Norweigen male population 1958-1997
#'
#' A dataset containing observed population numbers for Norweigen men.
#'
#' @format A data frame with 18 rows and 8 variables:
#' \describe{
#'   \item{rows}{Each row identifies a five year age interval, from 0-4 to 85+}
#'   \item{columns}{Each column identifies a five year period, from 1958-1962 to 1993-1997}
#' }
#' @examples
#' inpop1
#' cbind(inpop1, inpop2)
#' @family nordpred_example_data
"inpop1"



#' Predicted Norweigen male population 1998-2022
#'
#' A dataset containing predicted population numbers for Norweigen men.
#'
#' @format A data frame with 18 rows and 5 variables:
#' \describe{
#'   \item{rows}{Each row identifies a five year age interval, from 0-4 to 85+}
#'   \item{columns}{Each column identifies a five year period, from 1998-2002 to 2018-2022}
#' }
#' @examples
#' inpop2
#' cbind(inpop1, inpop2)
#' @family nordpred_example_data
"inpop2"
