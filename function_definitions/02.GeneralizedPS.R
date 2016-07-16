################################################################################
### Multinomial Logistic Regression to Add Generalized PS
## 
## Created on: 2015-03-16
## Author: Kazuki Yoshida
################################################################################

## Function to return a data frame with added PS
AddGPS <- function(data,
                   formula  = Tr ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                   psPrefix = "PS_",
                   family   = multinomial(parallel = FALSE)) {
    ## Fit multinomial logistic regression
    resVglm <- vglm(formula = formula,
                    data    = data,
                    family  = family)
    ## Calculate PS
    psData <- as.data.frame(predict(resVglm, type = "response"))
    names(psData) <- paste0(psPrefix, names(psData))
    ## Add to data
    cbind(data, psData)
}
