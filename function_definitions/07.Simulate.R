################################################################################
### Matching Weight and Algorithmic Matching Comparison
##
## Created on: 2015-03-17
## Author: Kazuki Yoshida
################################################################################


###
### Load necessary packages
################################################################################
library(magrittr)
library(VGAM)
library(survey)
library(doParallel)


###
### Load necessary functions
################################################################################

## Called from the directory above, thus, path needs to be the following
source("./function_definitions/02.GeneralizedPS.R")
source("./function_definitions/03.MatchingWeight.R")
source("./function_definitions/04.Rassen2013Frontend.R")
source("./function_definitions/05.BalanceCheck.R")
source("./function_definitions/06.FitModels.R")


###
### Function for one simulation iteration
################################################################################

###  Construct adjusted datasets given PS
ConstructDatasets <- function(data, twoPs, txVar, idVar,
                              txLevels, psPrefix,
                              minTrioCount, caliperFactor) {

    ## Calculate caliper as specified in Rassen 2013 p403
    caliper3Way <- CalcCaliper3Way(ps1 = data[,twoPs[1]],
                                   ps2 = data[,twoPs[2]],
                                   group = data[,txVar])

    ## Run Rassen 2013 matching method
    ## Loosen caliper by caliperFactor^i (where i = 0,1,2,...) until
    ## at least minTrioCount trios are created.
    ## Augmented data with matched id will return
    ## Requires minTrioCount <= nMinEach to be fail-safe
    ## Treatment groups have to be 1, 2, and 3
    out3WayM <- Match3WayAugmentData(data = data, caliper = caliper3Way,
                                     idVar = idVar, txVar = txVar,
                                     twoPs = twoPs,
                                     minTrioCount = minTrioCount,
                                     caliperFactor = caliperFactor)
    data <- out3WayM$outData

    ## Augment data with 3-group matching weights mw
    data <- MatchingWeightAugmentData(data = data,
                                      txVar = txVar,
                                      tx = txLevels,
                                      psPrefix = psPrefix)

    ## Generate adjusted datasets
    list(data      = data,
         matchData = subset(data, !is.na(set_num)),
         mwData    = svydesign(ids = ~ ID, weights = ~ mw,      data = data),
         iptwData  = svydesign(ids = ~ ID, weights = ~ st_iptw, data = data),
         out3WayM  = out3WayM)
}


###  Give true RRs based on counterfactual DRS and true PS
## The list has to have data, matchData, mwData, and iptwData
## The true counterfactual DRS has to be pY0, pY1, pY2.
TrueDRS <- function(lst){

    ## True RR for unmatched
    U <- c(mean(lst$data$pY0),
           mean(lst$data$pY1),
           mean(lst$data$pY2))

    ## True RR for matched
    M <- c(mean(lst$matchData$pY0),
           mean(lst$matchData$pY1),
           mean(lst$matchData$pY2))

    ## True RR for MW
    Mw <- c(as.numeric(svymean(~ pY0, lst$mwData)),
            as.numeric(svymean(~ pY1, lst$mwData)),
            as.numeric(svymean(~ pY2, lst$mwData)))

    ## True RR for IPTW
    Ip <- c(as.numeric(svymean(~ pY0, lst$iptwData)),
            as.numeric(svymean(~ pY1, lst$iptwData)),
            as.numeric(svymean(~ pY2, lst$iptwData)))

    trueDRS <- list(U = U, M = M, Mw = Mw, Ip = Ip)

    ## Name elements
    trueDRS <- lapply(trueDRS, function(vec) {
        names(vec) <- paste0("pY", c(0,1,2))
        vec
    })

    ## Obtain true RR
    trueCoef <- lapply(trueDRS, function(vec) {
        names(vec) <- NULL
        c(trueCoef.1v0 = log(vec[2] / vec[1]),
          trueCoef.2v0 = log(vec[3] / vec[1]),
          trueCoef.2v1 = log(vec[3] / vec[2]))
    })
    ##
    list(trueDRS,
         trueCoef)
}


###  Analyze one iteration
Simulate <- function(data, psFormula, lstParams) {

    ## Set Variable names
    idVar          <- "ID"
    outcomeVar     <- "Y"
    covPrefix      <- "X"
    covNames       <- paste0(covPrefix, 1:10)
    ## Binary covariate names (hard-coded)
    covBinary      <- paste0(covPrefix, c(4, 5, 10))
    ## Prefix for estimated propensity score variables
    psPrefix       <- "PS_"
    psTruePrefix   <- "pT"
    ## Caliper widening
    minTrioCount   <- 2
    caliperFactor  <- 1.5

    ## Derived variables
    txVar          <- as.character(psFormula[[2]])
    ## Effect modifier name (must be a binary variable)
    eModName       <- paste0(covPrefix, lstParams$eMod)
    stopifnot(eModName %in% covBinary)
    ## All treatment levels (assume numeric)
    txLevels       <- as.numeric(names(table(data[,txVar])))
    ## First two non-redundant PS to use
    twoPs          <- paste0(psPrefix, txLevels[1:2])
    twoPsTrue      <- paste0(psTruePrefix, txLevels[1:2])
    ## Outcome model formula
    outcomeFormula <- as.formula(sprintf("%s ~ factor(%s)", outcomeVar, txVar))


    ## Add ID done here to save storage for data
    data[,idVar] <- seq_len(nrow(data))

    ## Add GPS from multinomial logistic regression
    dataGps <- AddGPS(data, formula = psFormula, psPrefix = psPrefix)


    ## Generate adjusted datasets (U, M, Mw, Ip) using esimtated GPS
    lstDfs <- ConstructDatasets(data          = dataGps,
                                twoPs         = twoPs,
                                txVar         = txVar,
                                txLevels      = txLevels,
                                psPrefix      = psPrefix,
                                idVar         = idVar,
                                minTrioCount  = minTrioCount,
                                caliperFactor = caliperFactor)


    ## Generate adjusted datasets (U, M, Mw, Ip) using true GPS
    lstTrueDfs <- ConstructDatasets(data          = data,
                                    twoPs         = twoPsTrue,
                                    txVar         = txVar,
                                    txLevels      = txLevels,
                                    psPrefix      = psTruePrefix,
                                    idVar         = idVar,
                                    minTrioCount  = minTrioCount,
                                    caliperFactor = caliperFactor)
    ## Obtain true DRS for each dataset
    lstTrueDrs <- TrueDRS(lstTrueDfs)


    ## Common support of true GPS
    CS <- CommonSupport(data  = data,
                        txVar = txVar,
                        twoPs = twoPsTrue)
    ## Subset to common support by true GPS
    lstTrueDfsCs <- SubsetCommonSupport(CS        = CS,
                                        twoPs     = twoPsTrue,
                                        data      = lstTrueDfs$data,
                                        matchData = lstTrueDfs$matchData,
                                        mwData    = lstTrueDfs$mwData,
                                        iptwData  = lstTrueDfs$iptwData)
    ## Obtain true DRS for each dataset within common support
    lstTrueDrsCs <- TrueDRS(lstTrueDfsCs)
    ## Name each element ...Cs
    lstTrueDrsCs <- lapply(lstTrueDrsCs, function(lst) {
        ## Each list is a list of vectors
        lapply(lst, function(vec) {
            ## Change vector element names
            names(vec) <- gsub("trueCoef", "trueCoefCs", names(vec))
            vec
        })
    })



    ## Assessment begins here
    ## Assign datasets based on estimated PS for easier manipulation
    ## data is overwritten
    data      <- lstDfs$data
    matchData <- lstDfs$matchData
    mwData    <- lstDfs$mwData
    iptwData  <- lstDfs$iptwData
    out3WayM  <- lstDfs$out3WayM

    ## Sample sizes and caliper
    lstSizeCaliper <- list(U.n       = nrow(data),
                           ## Sample sizes
                           M.n       = nrow(matchData),
                           Mw.n      = sum(data$mw),
                           Ip.n      = sum(data$st_iptw),
                           ## Matching related information
                           ## Final caliper width
                           caliper   = out3WayM$caliper,
                           ## Power of caliper multiplication factor
                           caliperPw = out3WayM$caliperPower)

    ## Standardized mean differences
    lstSmds <- SmdAndMeanSmd(vars = covNames,
                             binaryVars = covBinary,
                             groupVar = txVar,
                             data = data, matchData = matchData,
                             mwData = mwData, iptwData = iptwData)

    ## Fit models and obtain a list of (coef, vcov) pairs
    lstModelOut <- FitModels(data = data, matchData = matchData,
                             mwData = mwData, iptwData = iptwData,
                             outcomeFormula = outcomeFormula)

    ## Summarize case (Y = 1) counts
    lstCaseCounts <- CaseCounts(data = data, matchData = matchData,
                                mwData = mwData, iptwData = iptwData,
                                outcomeVar = outcomeVar, txVar = txVar)

    ## Summarize X10 = 1 counts
    lstX10Counts <- CaseCounts(data = data, matchData = matchData,
                               mwData = mwData, iptwData = iptwData,
                               outcomeVar = "X10", txVar = txVar)
    vecX10Counts <- unlist(lstX10Counts)
    names(vecX10Counts) <- gsub("nCases", "X10", names(vecX10Counts))

    ## Summarize prevalence of effect modifier
    lstModifier <- PModifier(data = data, matchData = matchData,
                             mwData = mwData, iptwData = iptwData,
                             modifier = eModName, txVar = txVar)

    ## Summarize prevalence of effect modifier WITHIN COMMON SUPPORT
    lstModifierCs <- PModifierCs(data = data, matchData = matchData,
                                 mwData = mwData, iptwData = iptwData,
                                 modifier = eModName, txVar = txVar,
                                 twoPs = twoPs)

    ## Return results as an unlisted vector
    lstRes <- c(lstSizeCaliper,
                lstSmds,
                lstModelOut,
                lstTrueDrs,
                lstTrueDrsCs,
                lstCaseCounts,
                vecX10Counts,
                lstModifier,
                lstModifierCs)
    unlist(lstRes)
}



###  Variant to return adjusted datasets
SimulateReturnAdjData <- function(data, psFormula, lstParams) {

    ## Set Variable names
    idVar          <- "ID"
    outcomeVar     <- "Y"
    covPrefix      <- "X"
    covNames       <- paste0(covPrefix, 1:10)
    ## Binary covariate names (hard-coded)
    covBinary      <- paste0(covPrefix, c(4, 5, 10))
    ## Prefix for estimated propensity score variables
    psPrefix       <- "PS_"
    psTruePrefix   <- "pT"
    ## Caliper widening
    minTrioCount   <- 2
    caliperFactor  <- 1.5

    ## Derived variables
    txVar          <- as.character(psFormula[[2]])
    ## Effect modifier name (must be a binary variable)
    eModName       <- paste0(covPrefix, lstParams$eMod)
    stopifnot(eModName %in% covBinary)
    ## All treatment levels (assume numeric)
    txLevels       <- as.numeric(names(table(data[,txVar])))
    ## First two non-redundant PS to use
    twoPs          <- paste0(psPrefix, txLevels[1:2])
    twoPsTrue      <- paste0(psTruePrefix, txLevels[1:2])
    ## Outcome model formula
    outcomeFormula <- as.formula(sprintf("%s ~ factor(%s)", outcomeVar, txVar))


    ## Add ID done here to save storage for data
    data[,idVar] <- seq_len(nrow(data))

    ## Add GPS from multinomial logistic regression
    dataGps <- AddGPS(data, formula = psFormula, psPrefix = psPrefix)


    ## Generate adjusted datasets (U, M, Mw, Ip) using esimtated GPS
    lstDfs <- ConstructDatasets(data          = dataGps,
                                twoPs         = twoPs,
                                txVar         = txVar,
                                txLevels      = txLevels,
                                psPrefix      = psPrefix,
                                idVar         = idVar,
                                minTrioCount  = minTrioCount,
                                caliperFactor = caliperFactor)
    lstDfs
}


###
### Function for multiple iterations of one scenario
################################################################################

## Run multiple iterations of the same scenario given list of data frames
Iterate <- function(lstData, psFormula, lstParams) {

    simOut <- foreach(i = seq_along(lstData)) %dopar% {
        cat("### Working on dataset", i, "\n")
        Simulate(data = lstData[[i]],
                 psFormula = psFormula,
                 lstParams = lstParams)
    }

    as.data.frame(do.call(rbind, simOut))
}
