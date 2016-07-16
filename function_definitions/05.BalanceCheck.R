################################################################################
### Balance Checking Functions
##
## Created on: 2015-03-17
## Author: Kazuki Yoshida
################################################################################

### Load packages
library(magrittr)

### Custom tableone helper functions
source("./function_definitions/05.BalanceCheckHelpers.R")

### Function to assess SMDs for all datasets
SmdAndMeanSmd <- function(vars, binaryVars, groupVar,
                          data, matchData, mwData, iptwData) {

    ## For unmatched cohort
    SmdUnmatchRaw <- StdDiffs(data = data, vars = vars, groupVar = "Tr",
                                         binaryVars = binaryVars) %>% do.call(cbind, .)
    colnames(SmdUnmatchRaw) <- vars
    SmdUnmatch     <- colMeans(SmdUnmatchRaw)
    meanSmdUnmatch <- mean(SmdUnmatch)

    ## For matched cohort
    SmdMatchRaw <- StdDiffs(data = matchData, vars = vars, groupVar = "Tr",
                                       binaryVars = binaryVars) %>% do.call(cbind, .)
    colnames(SmdMatchRaw) <- vars
    SmdMatch     <- colMeans(SmdMatchRaw)
    meanSmdMatch <- mean(SmdMatch)

    ## For mw cohort
    SmdMwRaw <- svyStdDiffs(data = mwData, vars = vars, groupVar = "Tr",
                                       binaryVars = binaryVars) %>% do.call(cbind, .)
    colnames(SmdMwRaw) <- vars
    SmdMw     <- colMeans(SmdMwRaw)
    meanSmdMw <- mean(SmdMw)

    ## For iptw cohort
    SmdIpRaw <- svyStdDiffs(data = iptwData, vars = vars, groupVar = "Tr",
                                       binaryVars = binaryVars) %>% do.call(cbind, .)
    colnames(SmdIpRaw) <- vars
    SmdIp     <- colMeans(SmdIpRaw)
    meanSmdIp <- mean(SmdIp)

    ## Return as a vector named appropriately
    list(U.smd      = SmdUnmatch,
         U.meanSmd  = meanSmdUnmatch,

         M.smd      = SmdMatch,
         M.meanSmd  = meanSmdMatch,

         Mw.smd     = SmdMw,
         Mw.meanSmd = meanSmdMw,

         Ip.smd     = SmdIp,
         Ip.meanSmd = meanSmdIp)
}
