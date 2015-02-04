## ------------------------------------------------------------------------
library(nordpred)

## ------------------------------------------------------------------------
knitr::kable(inpop1)

## ------------------------------------------------------------------------
knitr::kable(inpop2)

## ------------------------------------------------------------------------
knitr::kable(inpop2)

## ------------------------------------------------------------------------
inpop <- cbind(inpop1, inpop2)

## ------------------------------------------------------------------------
est <- nordpred.estimate(cases=indata, pyr = inpop, noperiod = 4, startestage = 5)

## ------------------------------------------------------------------------
res <- nordpred.prediction(est, startuseage = 6, cuttrend = c(0,.25,.5,.75,.75), recent=TRUE)

## ------------------------------------------------------------------------
res <- nordpred(cases = indata, pyr = inpop, startestage = 5, startuseage = 6, 
                noperiods = 4, cuttrend = c(0, .25, .5, .75, .75))

## ------------------------------------------------------------------------
res <- nordpred(indata, inpop, startestage = 5, startuseage = 6, noperiods = 4:6,
                cuttrend=c(0, .25, .5, .75, .75))


## ------------------------------------------------------------------------
est2 <- nordpred.estimate(indata, inpop, 4, 5, linkfunc = "poisson")
res2 <- nordpred.prediction(est2, startuseage = 6, cuttrend = c(0, .25, .5, .75, .75), recent = TRUE)

## ------------------------------------------------------------------------
print(res)
nordpred.getpred(res)
summary(res, printpred = FALSE)

## ------------------------------------------------------------------------
## World population standard
wstand <- c(0.12, 0.1, 0.09, 0.09, 0.08, 0.08, 0.06, 0.06, 0.06, 0.06,0.05, 
            0.04, 0.04, 0.03, 0.02, 0.01, 0.005, 0.005)
            
round(nordpred.getpred(res, incidence = TRUE, standpop = NULL), 2)
round(nordpred.getpred(res, incidence = TRUE, standpop = wstand), 2)


## ------------------------------------------------------------------------
plot(res, standpop = wstand)

## ------------------------------------------------------------------------
plot(res2, standpop = wstand)
plot(res, new = FALSE, lty = c(1, 2), standpop = wstand)

## ------------------------------------------------------------------------
plot(nordpred.prediction(est, startuseage = 6, cuttrend = c(0, 0, 0, 0, 0), 
    recent = FALSE), standpop = wstand, new = TRUE)
plot(nordpred.prediction(est, startuseage = 6, cuttrend = c(1, 1, 1, 1, 1), 
    recent = FALSE), standpop = wstand, new = FALSE, lty = c(1, 2))
plot(nordpred.prediction(est, startuseage = 6, cuttrend = c(0, .25, .5, .75, .75), 
    recent = FALSE), standpop = wstand, new = FALSE, lty = c(1, 4))

