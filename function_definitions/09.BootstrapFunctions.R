################################################################################
### Functions for bootstrapping
##
## Created on: 2016-01-19
## Author: Kazuki Yoshida
################################################################################


library(boot)

source("./function_definitions/07.Simulate.R")


### Bootstrap function constructor
ConstructBootFun <- function(psFormula, lstParams) {

    ## Return a function appropriate for boot()
    function(data, ii) {
        ## Resample
        data <- data[ii,]

        ## Set Variable names
        idVar          <- "ID"
        outcomeVar     <- "Y"
        covPrefix      <- "X"
        covNames       <- paste0(covPrefix, 1:10)
        ## Prefix for estimated propensity score variables
        psPrefix       <- "PS_"
        ## Derived variables
        txVar          <- as.character(psFormula[[2]])
        ## All treatment levels (assume numeric)
        txLevels       <- as.numeric(names(table(data[,txVar])))
        ## First two non-redundant PS to use
        twoPs          <- paste0(psPrefix, txLevels[1:2])

        ## Add ID done here to save storage for data
        data[,idVar] <- seq_len(nrow(data))

        ## Add GPS from multinomial logistic regression
        dataGps <- AddGPS(data, formula = psFormula, psPrefix = psPrefix)

        ## Add weights
        dataGps <- MatchingWeightAugmentData(data = dataGps,
                                             txVar = txVar,
                                             tx = txLevels,
                                             psPrefix = psPrefix)

        ## Create a weighted data object
        mwData <- svydesign(ids = ~ ID, weights = ~ mw, data = dataGps)

        ## Fit outcome model
        outcomeFormula <- as.formula(sprintf("%s ~ factor(%s)", outcomeVar, txVar))
        modelMw <- svyglm(formula = outcomeFormula, design = mwData, family = quasipoisson())

        ## Coefs for intercept, 1v0, 2v0
        coefs <- coef(modelMw)
        ## Add 2v1
        coefs[4] <- coefs[3] - coefs[2]
        ## Drop intercept
        coefs <- coefs[-1]
        ## Rename
        names(coefs) <- c("1v0","2v0","2v1")
        coefs
    }
}
