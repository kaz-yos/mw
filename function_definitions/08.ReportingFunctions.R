################################################################################
### Functions used for reporting
##
## Created on: 2015-07-15
## Author: Kazuki Yoshida
################################################################################


library(magrittr)

###
### Summary result functions
###################################################################################

### String converter with rounding
r <- function(x, digits = 3) {
    sprintf(fmt = paste0("%.", digits,"f"), x)
}

### Summarize coverage
## This takes a p vector and two nxp matrices
Coverage <- function(trueVec, lowerMat, upperMat) {
    stopifnot(ncol(lowerMat) == length(trueVec) & ncol(upperMat) == length(trueVec))

    ## TRUE if lower < true
    grThanLower <- sweep(x = lowerMat, MARGIN = 2, STATS = trueVec,
                         FUN = "<") # trueVec on RHS
    ## TRUE if upper > true
    lsThanUpper <- sweep(x = upperMat, MARGIN = 2, STATS = trueVec,
                         FUN = ">") # trueVec on RHS
    ## Proportion of meeting both
    coverage <- colMeans(grThanLower * lsThanUpper)
    names(coverage) <- gsub("lower", "cvr", names(coverage))
    ## Return a vector
    coverage
}


###
### Aggregate results from individual summary function
################################################################################

### Function to pick up names matching pat
PickMatchNames <- function(chrVec, pat) {
    Filter(f = function(elt) {grepl(pat, elt)}, x = chrVec)
}

###
### Report one scenario
Report <- function(scenarioCount, lstParams, dfOut) {

    ## Sort variable names by methods
    varNamesUnsort <- names(dfOut)
    varNames <- c(PickMatchNames(chrVec = varNamesUnsort, pat = "^U\\."),
                  PickMatchNames(chrVec = varNamesUnsort, pat = "^M\\."),
                  PickMatchNames(chrVec = varNamesUnsort, pat = "^Mw\\."),
                  PickMatchNames(chrVec = varNamesUnsort, pat = "^Ip\\."),
                  PickMatchNames(chrVec = varNamesUnsort, pat = "caliper"))

    ## Pick variable names that just need averaging
    varsToAvg <- c(PickMatchNames(chrVec = varNames, pat = "\\.n$"),
                   ## Sample size in common support
                   PickMatchNames(chrVec = varNames, pat = "\\.nCs$"),
                   ## SMDs
                   PickMatchNames(chrVec = varNames, pat = "\\.smd\\."),
                   PickMatchNames(chrVec = varNames, pat = "\\.meanSmd"),
                   ## Point estimates (averaging on log scale)
                   PickMatchNames(chrVec = varNames, pat = "\\.coef"),
                   ## True coefficients (averaging on log scale)
                   PickMatchNames(chrVec = varNames, pat = "\\.trueCoef"),
                   ## Estimated variances
                   PickMatchNames(chrVec = varNames, pat = "\\.vars"),
                   ## Caliper related
                   PickMatchNames(chrVec = varNames, pat = "caliper"),
                   ## Modifier prevalence (including pModCs in common support)
                   PickMatchNames(chrVec = varNames, pat = "pMod"),
                   ## Case count
                   PickMatchNames(chrVec = varNames, pat = "nCases"),
                   ## (X10 = 1) count
                   PickMatchNames(chrVec = varNames, pat = "[UMMwIp].X10"))

    ## Do not remove NaN! It fails when it should.
    means <- colMeans(dfOut[varsToAvg], na.rm = FALSE)


### True variance calculation using coefficients
    trueVars <- apply(dfOut[,PickMatchNames(chrVec = varNames, pat = "coef")],
                      MARGIN = 2, FUN = var)
    names(trueVars) <- gsub("coef", "trueV", names(trueVars))


### Rejection rate (Power in non-null betaT scenario; alpha in null betaT scenario)
    ## z can be negative or positive
    rejRate <- colMeans(abs(dfOut[,PickMatchNames(chrVec = varNames, pat = "\\.z\\.")]) > 1.95996)
    names(rejRate) <- gsub("z", "pRej", names(rejRate))


### True marginal RR calculation based on true log RR
    ## Average iteration-specific true log-RR, and the exponentiate
    ## Estimated coefficients are also averaged, and then exponentiated
    expectMargRR <- exp(means[PickMatchNames(chrVec = varNames, pat = "trueCoef\\.")])
    names(expectMargRR) <- gsub("trueCoef", "trueRR", names(expectMargRR))
    ## Estimate expected marginal RR within common support
    expectMargRRCs <- exp(means[PickMatchNames(chrVec = varNames, pat = "trueCoefCs\\.")])
    names(expectMargRRCs) <- gsub("trueCoef", "trueRR", names(expectMargRRCs))


### Bias calculation
    ## Mean estimated coefficients exponentiated (exp(mean(coef)))
    ## Averaging came first
    vecEstMarRR <- exp(means[PickMatchNames(chrVec = varNames, pat = "coef")])
    names(vecEstMarRR) <- gsub("coef", "RR", names(vecEstMarRR))
    ## (Mean estimated marginal RR) / (expected marginal RR)
    estRRDivExpMargRR <- vecEstMarRR / expectMargRR
    names(estRRDivExpMargRR) <- gsub("RR", "biasRR", names(estRRDivExpMargRR))
    ## Using true RR within common support
    estRRDivExpMargRRCs <- vecEstMarRR / expectMargRRCs
    names(estRRDivExpMargRRCs) <- gsub("RR", "biasRRCs", names(estRRDivExpMargRRCs))


### Mean squared error calculation
    ## Calculate for coefficient (log RR) averaged across iterations
    ## MSE = var + bias^2
    mse   <- trueVars + (means[PickMatchNames(chrVec = varNames, pat = "coef")] - log(expectMargRR))^2
    names(mse) <- gsub("trueV", "mse", names(mse))
    ## Based on common support true marginal log RR
    mseCs <- trueVars + (means[PickMatchNames(chrVec = varNames, pat = "coef")] - log(expectMargRRCs))^2
    names(mseCs) <- gsub("trueV", "mseCs", names(mseCs))


### Coverage of the expected true log RR
    coverage <- Coverage(true  = log(expectMargRR),
                         lower = dfOut[,PickMatchNames(chrVec = varNames, pat = "lower")],
                         upper = dfOut[,PickMatchNames(chrVec = varNames, pat = "upper")])
    ## Using true RR within common support
    coverageCs <- Coverage(true  = log(expectMargRRCs),
                           lower = dfOut[,PickMatchNames(chrVec = varNames, pat = "lower")],
                           upper = dfOut[,PickMatchNames(chrVec = varNames, pat = "upper")])
    names(coverageCs) <- gsub("cvr", "cvrCs", names(coverageCs))


### Minimum counts
    minCaseCount <- apply(dfOut[,PickMatchNames(chrVec = varNames, pat = "nCases")],
                          MARGIN = 2, FUN = min)
    names(minCaseCount) <- gsub("nCases", "minNCases", names(minCaseCount))

    minX10Count <- apply(dfOut[,PickMatchNames(chrVec = varNames, pat = "[UMMwIp].X10")],
                         MARGIN = 2, FUN = min)
    names(minX10Count) <- gsub("X10", "minX10", names(minX10Count))

    ## Return the entire thing as a vector
    c(scenario = scenarioCount,
      unlist(lstParams),
      ## Results that only required averaging
      means,
      ## True variance based on simulation
      trueVars,
      ## Rejection rate (power or alpha)
      rejRate,
      ## Expected true marginal RR based on true counterfactual DRS
      expectMargRR,
      ## Expected true marginal RR within common support
      expectMargRRCs,
      ## Bias in marginal RR: mean(estimated RR) / expected RR
      estRRDivExpMargRR,
      ## Bias in marginal RR using true marginal RR within common support
      estRRDivExpMargRRCs,
      ## MSE of coefficients
      mse,
      ## MSE of coefficients within common support
      mseCs,
      ## Coverage probability
      coverage,
      ## Coverage probability using true marginal RR within common support
      coverageCs,
      ## Minimum case count over iterations
      minCaseCount,
      ## Minimum X10 count over iterations
      minX10Count)
}


### Aggregate report
Gather <- function(data, indexVars, gatherVarPat) {

    gatherVars <- PickMatchNames(names(data), gatherVarPat)

    out <- gather_(data = data[, c(indexVars, gatherVars)],
                   key_col = "key", value_col = "value",
                   gather_cols = gatherVars)

    out$method   <- factor(gsub("\\..*", "", as.character(out$key)), levels = c("U","M","Mw","Ip"))
    out$contrast <- factor(gsub(".*\\.", "", as.character(out$key)), levels = c("1v0","2v0","2v1"))

    out
}
