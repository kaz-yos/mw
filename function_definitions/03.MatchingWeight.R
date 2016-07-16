################################################################################
### Create matching weight given generalized propensity scores
##
## Created on: 2015-07-22
## Author: Kazuki Yoshida
################################################################################

### Function to augment the data frame with matching weight mw
MatchingWeightAugmentData <- function(data, txVar = "Tr", tx = c(0,1,2), psPrefix = "PS_") {

    ## Treatment indicator data frame (any number of groups allowed)
    dfAssign <- as.data.frame(lapply(tx, function(tx_k) {
        as.numeric(data[txVar] == tx_k)
    }))
    colnames(dfAssign) <- paste0(txVar, tx)

    ## Name of PS variables
    psVars <- paste0(psPrefix, tx)

    ## Pick denominator (PS for assigned treatment)
    data$PSassign <- rowSums(data[psVars] * dfAssign)

    ## Calculate iptw
    data$iptw <- exp(- log(data$PSassign))

    ## Calculate marginal prevalence of assigned treatment
    ## This assumes tx argument is correctly sorted.
    txProp <- prop.table(table(data[,txVar]))
    txPropAssign <- as.numeric(as.matrix(dfAssign) %*% txProp)

    ## Calculate stabilized iptw
    data$st_iptw <- exp(log(txPropAssign) - log(data$PSassign))

    ## Pick numerator for MW (smallest of all PS)
    data$PSmin <- do.call(pmin, data[psVars])

    ## Calculate matching weight
    data$mw <- exp(log(data$PSmin) - log(data$PSassign))

    ## Return the whole data
    data
}
