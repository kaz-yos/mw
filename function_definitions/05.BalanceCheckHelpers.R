################################################################################
### Balance check helpers (custom functions extracted from tableone)
##
## Created on: 2016-01-18
## Author: Kazuki Yoshida
################################################################################


###
### Helpers common to both unweighted and weighed functions
################################################################################


###
### Functions for unweighted data only
################################################################################

### Check strata for NA-only strata
CheckNaOnlyStrata <- function(variable, group) {

    unlist(lapply(split(variable, group),
                  function(var) {
                      ## TRUE if only NA's within stratum
                      all(is.na(var))
                  }))
}


### Continuous/binary standardized mean differences
## Expects continuous or 0,1 binary variable
StdDiff <- function(variable, group, binary = FALSE, na.rm = TRUE) {

    ## Check strata for all NA strata
    logiAllNaStrata <- CheckNaOnlyStrata(variable, group)
    ## If ANY strata have only NA's do not remove NA's
    if (any(logiAllNaStrata)) {
        warning("Variable has only NA's in at least one stratum. na.rm turned off.")
        na.rm = FALSE
    }

    ## Proportion of 1 is the mean of variable
    means <- tapply(variable, group, mean, na.rm = na.rm)

    ## Variance is p(1-p)
    if (binary) {
        vars  <- means * (1 - means)
    } else {
        vars  <- tapply(variable, group, var, na.rm = na.rm)
    }

    ## Outer to obtain all pairwise differences
    meanDiffs  <- outer(X = means, Y = means, FUN = "-")
    ## Outer to obtain all pairwise variance mean
    varMeans   <- outer(X = vars, Y = vars, FUN = "+") / 2

    out <- meanDiffs / sqrt(varMeans)

    ## If mean difference is zero and variance is zero, these are the same constants.
    ## In this case, differences should be defined zero.
    out[is.na(out) & (meanDiffs == 0) & (varMeans == 0)] <- 0

    ## This lower.tri() approach is actually giving 2vs1, 3vs1, etc
    ## opposite of stated 1vs2, 1vs3. Only correct if abs() is used.
    abs(out[lower.tri(out)])
}


### Standardized mean differences for multiple variables
## Continuous or binary only
StdDiffs <- function(data, vars, groupVar, binaryVars) {

    lapply(vars, function(var) {

        StdDiff(variable = data[,var], group = data[,groupVar],
                binary = (var %in% binaryVars))
    })
}


###
### Functions for weighted data only
################################################################################

### Check strata for NA-only strata
svyCheckNaOnlyStrata <- function(varName, groupName, design) {

    unlist(lapply(split(design$variables[,varName],
                        design$variables[,groupName]),
                  function(var) {
                      ## TRUE if only NA's within stratum
                      all(is.na(var))
                  }))
}

### Continuous/binary standardized mean differences
## Expects continuous or 0,1 binary variable
svyStdDiff <- function(varName, groupName, design, binary = FALSE, na.rm = TRUE) {

    ## Check strata for all NA strata
    logiAllNaStrata <- svyCheckNaOnlyStrata(varName, groupName, design)
    ## If ANY strata have only NA's do not remove NA's
    if (any(logiAllNaStrata)) {
        warning(varName, " has only NA's in at least one stratum. na.rm turned off.")
        na.rm = FALSE
    }

    varFormula   <- as.formula(paste("~", varName))
    groupFormula <- as.formula(paste("~", groupName))

    means <- svyby(formula = varFormula, by = groupFormula,
                   FUN = svymean, design = design, na.rm = na.rm)[,2]

    if (binary) {
        vars  <- means * (1 - means)
    } else {
        vars  <- svyby(formula = varFormula, by = groupFormula,
                       FUN = svyvar,  design = design, na.rm = na.rm)[,2]
    }

    ## Outer to obtain all pairwise differences
    meanDiffs <- outer(X = means, Y = means, FUN = "-")
    ## Outer to obtain all pairwise variance mean
    varMeans  <- outer(X = vars, Y = vars, FUN = "+") / 2

    out <- meanDiffs / sqrt(varMeans)

    ## If mean difference is zero and variance is zero, these are the same constants.
    ## In this case, differences should be defined zero.
    out[is.na(out) & (meanDiffs == 0) & (varMeans == 0)] <- 0

    ## This lower.tri() approach is actually giving 2vs1, 3vs1, etc
    ## opposite of stated 1vs2, 1vs3. Only correct if abs() is used.
    abs(out[lower.tri(out)])
}


### Standardized mean differences for multiple variables
## Continuous or binary only
svyStdDiffs <- function(data, vars, groupVar, binaryVars) {

    lapply(vars, function(var) {

        svyStdDiff(varName = var, groupName = groupVar, design = data,
                   binary = (var %in% binaryVars))
    })
}
