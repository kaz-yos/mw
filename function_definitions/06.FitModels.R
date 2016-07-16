################################################################################
### Model fitting functions
##
## Created on: 2015-12-24
## Author: Kazuki Yoshida
################################################################################

## Load sandwich package for robust sandwich covariance matrix estimators
library(sandwich)
library(dplyr)


### Add confidence interval
AddCiAndZ <- function(ciQuantile = 1.95996, lstCoefVar) {

    coef  <- lstCoefVar$coef
    vars  <- lstCoefVar$vars
    se    <- sqrt(vars)
    lower <- coef - ciQuantile * se
    upper <- coef + ciQuantile * se
    z     <- coef / se

    ## Return minimum results
    list(coef = coef, vars = vars, lower = lower, upper = upper, z = z)
}


### Create all contrasts and then add CI and Z value
ThreeContrasts <- function(lstCoefVcov) {
    coef <- lstCoefVcov$coef
    ## Add third contrasts
    coef[4] <- coef[3] - coef[2]
    ## Name contrasts
    names(coef) <- c("Int","1v0","2v0","2v1")

    ## Variance for third contrast
    vcov <- lstCoefVcov$vcov
    vars <- diag(vcov)
    vars[4] <- vcov[2,2] + vcov[3,3] - 2 * vcov[2,3]
    names(vars) <- names(coef)

    ## Drop intercept at position 1
    lstCoefVar <- list(coef = coef[-1], vars = vars[-1])

    ## Return after adding CI and Z value
    AddCiAndZ(ciQuantile = 1.95996, lstCoefVar)
}


### Fit models and obtain correct estimates
FitModels <- function(data, matchData, mwData, iptwData,
                      outcomeFormula = Y ~ factor) {
    ## Unadjusted
    modelUnmatch <- glm(formula = outcomeFormula, data = data,       family = quasipoisson())
    ## Three-way Matching
    modelMatch   <- glm(formula = outcomeFormula, data = matchData,  family = quasipoisson())
    ## MW
    modelMw   <- svyglm(formula = outcomeFormula, design = mwData,   family = quasipoisson())
    ## IPTW
    modelIptw <- svyglm(formula = outcomeFormula, design = iptwData, family = quasipoisson())

    ## Obtain appropriate robust sandwich variance-covariance matrices
    ## svyglm already uses sandwich estimator, just use vcov
    ## https://stat.ethz.ch/pipermail/r-help/2008-December/183272.html

    ## Return pairs of coef and vcov
    lstOut <- list(U  = list(coef = coef(modelUnmatch), vcov = sandwich(modelUnmatch)),
                   M  = list(coef = coef(modelMatch),   vcov = sandwich(modelMatch)),
                   Mw = list(coef = coef(modelMw),      vcov = vcov(modelMw)),
                   Ip = list(coef = coef(modelIptw),    vcov = vcov(modelIptw)))

    ## Return as one list containing coef, vars, CI, and Z for each model
    sapply(lstOut, ThreeContrasts,
           simplify = FALSE)

}


### Obtain prevalence of effect modifier
PModifier <- function(data, matchData, mwData, iptwData,
                      modifier = "X5", txVar = "Tr") {

    form1  <- as.formula(sprintf(" ~ %s", modifier))
    formBy <- as.formula(sprintf(" ~ %s", txVar))

    U.pMod <- c(pMod = mean(data[,modifier]), pMod = tapply(data[,modifier], data[,txVar], FUN = mean))
    M.pMod <- c(pMod = mean(matchData[,modifier]), pMod = tapply(matchData[,modifier], matchData[,txVar], FUN = mean))

    Mw.pMod <- c(as.numeric(svymean(form1, design = mwData)),
                 svyby(form1, by = ~ Tr, design = mwData, FUN = svymean)[,modifier])
    names(Mw.pMod) <- names(U.pMod)
    Ip.pMod <- c(as.numeric(svymean(form1, design = iptwData)),
                 svyby(form1, by = ~ Tr, design = iptwData, FUN = svymean)[,modifier])
    names(Ip.pMod) <- names(U.pMod)

    list(U  = U.pMod,
         M  = M.pMod,
         Mw = Mw.pMod,
         Ip = Ip.pMod)
}

### Obtain ranges of common support
CommonSupport <- function(data, txVar = "Tr", twoPs = c("PS_0","PS_1")) {
    ## Within group extremes
    mins <- data[,c(txVar, twoPs)] %>%
        group_by(Tr) %>%
        summarize_each(funs(min))

    maxs <- data[,c(txVar, twoPs)] %>%
        group_by(Tr) %>%
        summarize_each(funs(max))

    ## Max of mins
    maxOfMins<- apply(mins[,-1], MARGIN = 2, FUN = max)
    ## Min of maxs
    minOfMaxs <- apply(maxs[,-1], MARGIN = 2, FUN = min)
    ##
    list(maxOfMins = maxOfMins, minOfMaxs = minOfMaxs)
}

### Subset datasets given common support information
SubsetCommonSupport <- function(CS, twoPs,
                                data, matchData, mwData, iptwData) {
    ## Subset to common support
    ## original cohort
    indWithinCs  <-
        CS$maxOfMins[1] <= data[,twoPs[1]] & data[,twoPs[1]] <= CS$minOfMaxs[1] &
        CS$maxOfMins[2] <= data[,twoPs[2]] & data[,twoPs[2]] <= CS$minOfMaxs[2]
    dataSub      <- data[indWithinCs, ]
    ## Matched cohort (This is trimming after matching)
    indWithinCsM <-
        CS$maxOfMins[1] <= matchData[,twoPs[1]] & matchData[,twoPs[1]] <= CS$minOfMaxs[1] &
        CS$maxOfMins[2] <= matchData[,twoPs[2]] & matchData[,twoPs[2]] <= CS$minOfMaxs[2]
    matchDataSub <- matchData[indWithinCsM, ]
    ## MW: All patients are involved, the same indices as the unmatched can be used.
    mwDataSub    <- subset(mwData, indWithinCs)
    ## IPTW: All patients are involved, the same indices as the unmatched can be used.
    iptwDataSub  <- subset(iptwData, indWithinCs)

    list(data      = dataSub,
         matchData = matchDataSub,
         mwData    = mwDataSub,
         iptwData  = iptwDataSub)
}

### Obtain prevalence of effect modifier within the common support
PModifierCs <- function(data, matchData, mwData, iptwData,
                        modifier = "X5",
                        txVar = "Tr", twoPs = c("PS_0","PS_1")) {

    ## Obtain the common support
    CS <- CommonSupport(data = data, txVar = txVar, twoPs = twoPs)

    ## Subset to common support
    ## original cohort
    indWithinCs  <-
        CS$maxOfMins[1] <= data[,twoPs[1]] & data[,twoPs[1]] <= CS$minOfMaxs[1] &
        CS$maxOfMins[2] <= data[,twoPs[2]] & data[,twoPs[2]] <= CS$minOfMaxs[2]
    dataSub      <- data[indWithinCs, ]
    ## Matched cohort
    indWithinCsM <-
        CS$maxOfMins[1] <= matchData[,twoPs[1]] & matchData[,twoPs[1]] <= CS$minOfMaxs[1] &
        CS$maxOfMins[2] <= matchData[,twoPs[2]] & matchData[,twoPs[2]] <= CS$minOfMaxs[2]
    matchDataSub <- matchData[indWithinCsM, ]
    ## MW: All patients are involved, the same indices as the unmatched can be used.
    mwDataSub    <- subset(mwData, indWithinCs)
    ## IPTW: All patients are involved, the same indices as the unmatched can be used.
    iptwDataSub  <- subset(iptwData, indWithinCs)

    ## Obtain sample size with in the common support
    lstN <- list(U  = sum(indWithinCs),
                 M  = sum(indWithinCsM),
                 ## allprob has inverse of weights as a data frame
                 Mw = sum(1 / mwDataSub$allprob[,1]),
                 Ip = sum(1 / iptwDataSub$allprob[,1]))

    ## Obtain marginal prevalence of modifier in the common support
    ## Just pass the subsetted datasets to PModifier()
    lstPMCs <- PModifier(data      = dataSub,
                         matchData = matchDataSub,
                         mwData    = mwDataSub,
                         iptwData  = iptwDataSub,
                         modifier  = modifier,
                         txVar     = txVar)

    ## Make sure they have the same lengths
    stopifnot(length(lstPMCs) == length(lstN))

    ## Rename the vectors in the output list appropriately
    ## Add n information
    lstOut <- sapply(names(lstPMCs), function(name) {
        ## Extract data
        out <- lstPMCs[[name]]
        names(out) <- gsub("pMod", "pModCs", names(out))
        ## Add n in common support
        c(out, nCs = lstN[[name]])
    }, simplify = FALSE)
    lstOut
}


### Obtain case counts overall and within treatment category
CaseCounts <- function(data, matchData, mwData, iptwData,
                       outcomeVar = "Y", txVar = "Tr"){

    formOut <- as.formula(sprintf(" ~ %s", outcomeVar))
    formBy  <- as.formula(sprintf(" ~ %s", txVar))

    U.nCases <- c(nCases = sum(data[,outcomeVar]),
                  nCases = tapply(data[,outcomeVar], data[,txVar], FUN = sum))

    M.nCases <- c(nCases = sum(matchData[,outcomeVar]),
                  nCases = tapply(matchData[,outcomeVar], matchData[,txVar], FUN = sum))

    Mw.nCases <- c(as.numeric(svytotal(formOut, design = mwData)),
                   svyby(formOut, by = ~ Tr, design = mwData, FUN = svytotal)[,outcomeVar])
    names(Mw.nCases) <- names(U.nCases)

    Ip.nCases <- c(as.numeric(svytotal(formOut, design = iptwData)),
                   svyby(formOut, by = ~ Tr, design = iptwData, FUN = svytotal)[,outcomeVar])
    names(Ip.nCases) <- names(U.nCases)

    list(U  = U.nCases,
         M  = M.nCases,
         Mw = Mw.nCases,
         Ip = Ip.nCases)
}
