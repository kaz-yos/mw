################################################################################
### Functions for Franklin 2014 data generation
##
## Created on: 2015-12-24
## Author: Kazuki Yoshida
################################################################################

## Franklin 2014 method for the covariate matrix X is used
CreateX <- function(N) {

    ## Continuous
    X1     <- rnorm(n = N, mean = 0, sd = 1)
    X2     <- rlnorm(n = N, meanlog = 0, sdlog = 0.5)
    X3     <- rnorm(n = N, mean = 0, sd = 10)

    ## Binary
    ##  Dependent on X1; mean(pX4) = 0.5
    oddsX4 <- exp(2 * X1)
    pX4    <- oddsX4 / (1 + oddsX4)
    X4     <- rbinom(n = N, size = 1, prob = pX4)
    ##  Independent
    X5     <- rbinom(n = N, size = 1, prob = 0.2)

    ## Multinomial
    X6     <- as.numeric(c(1:5) %*% rmultinom(n = N, size = 1, prob = c(0.5,0.3,0.1,0.05,0.05)))

    ## Induced variables
    X7     <- sin(X1)
    X8     <- X2^2
    X9     <- X3 * X4
    X10    <- X4 * X5

    ## Output a matrix
    cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
}

## mu vector for E[Xi]
## Numerically estimate for X4
## Use lnN(mu,sigma^2) for X2 and X8
## Use independence for X9 and X10
mu <- c(X1 = 0, X2 = exp(0.5^2/2), X3 = 0, X4 = 0.5, X5 = 0.2,
        X6 = 1.85, X7 = 0, X8 = exp(2*(0.5)^2), X9 = 0*0.5, X10 = 0.5*0.2)

## Generate N multinomial variables given NxK probability matrix
rmultinomMat <- function(pMat) {
    ## Convert to a KxN data frame
    pLst <- as.data.frame(t(pMat))

    ## Loop over individuals
    lstT <- lapply(pLst, function(pVec) {
        c(0,1,2) %*% rmultinom(n = 1, size = 1, prob = pVec)
    })

    ## Vector of individual treatment assignment
    unlist(lstT)
}

## Treatment assignment mechanism given a multinomial logistic model
AssignT <- function(alpha10, alpha20, alpha1X, alpha2X, X) {
    ## Log probability ratio: log(P(T=1)/P(T=0))
    logPr1v0 <- alpha10 +  X %*% alpha1X
    ## Log probability ratio: log(P(T=2)/P(T=0))
    logPr2v0 <- alpha20 +  X %*% alpha2X

    ## Probability ratio: P(T=1)/P(T=0)
    pr1v0    <- exp(logPr1v0)
    ## Probability ratio: P(T=2)/P(T=0)
    pr2v0    <- exp(logPr2v0)

    ## Normalizing constant
    denom    <- 1 + pr1v0 + pr2v0

    ## Individual treatment probabilities
    p0       <- 1     / denom
    p1       <- pr1v0 / denom
    p2       <- pr2v0 / denom
    ## NxK matrix (this configuration for df conversion)
    pMat     <- cbind(p0, p1, p2)

    ## Sample treatment assignments
    ## This requires N invocations and is inefficient.
    Tr <- rmultinomMat(pMat)

    ## Return various components together
    data.frame(pT0 = p0,                 # True propensity score for Tr = 0
               pT1 = p1,                 # True propensity score for Tr = 1
               pT2 = p2,                 # True propensity score for Tr = 2
               Tr  = Tr)                 # Trinomial treatment assignment
}

##
TransformLp <- function(lp) {
    ## Log-probability model
    ## Use lp as log probability and do exp(.) transformaation (log-linear model)
    expLp <- exp(lp)
    ## Trim probability greater than 1
    pY <- pmin(expLp, 1)

    ## Logistic model
    ## Use lp as log odds and do expit(.) transfromation (logistic model)
    ## oddsY <- exp(lp)
    ## pY <- oddsY / (1 + oddsY)

    ## pY in [0,1]
    stopifnot(all(0 <= pY & pY <=1))

    ## Return with exponentiated lp
    list(expLp = expLp, pY = pY)
}

## Outcome assignment mechanism given a log-linear model
## Assumes Tr in {0,1,2}
AssignY <- function(beta0, betaX, betaT, betaTX, eMod, X, Tr) {

    ## Detect unexpected treatment values
    stopifnot(all(Tr %in% c(0,1,2)))
    ## Detect non-binary effect modifiers
    stopifnot(all(X[,eMod] %in% c(0,1)))

    ## Counterfactual linear predictors
    ## sum(c(0,0) * betaT) and sum(c(0,0) * betaTX) give appropriate scalars.
    lp0  <- beta0 + (X %*% betaX) + sum(c(0,0) * betaT) + X[,eMod] * sum(c(0,0) * betaTX)
    lp1  <- beta0 + (X %*% betaX) + sum(c(1,0) * betaT) + X[,eMod] * sum(c(1,0) * betaTX)
    lp2  <- beta0 + (X %*% betaX) + sum(c(0,1) * betaT) + X[,eMod] * sum(c(0,1) * betaTX)

    ## Treatment dummy variables
    TMatAll <- cbind(as.numeric(Tr == 0),
                     as.numeric(Tr == 1),
                     as.numeric(Tr == 2))

    ## One treatment for all inviduals
    stopifnot(all(rowSums(TMatAll) == 1))

    ## Assing linear predictor based on treatment assignment
    lp <- rowSums(cbind(lp0, lp1, lp2) * TMatAll)

    ## Transform to P(Y=1)
    transformLp <- TransformLp(lp)
    ## Counterfactual versions
    pY0  <- TransformLp(lp0)$pY
    pY1  <- TransformLp(lp1)$pY
    pY2  <- TransformLp(lp2)$pY

    ## Sample binary Y
    Y    <- rbinom(n = length(transformLp$pY), size = 1, prob = transformLp$pY)

    ## Return various components together
    data.frame(lp    = lp,                  # Linear predictor
               expLp = transformLp$expLp,   # Exponential of lp
               pY    = transformLp$pY,      # True DRS
               pY0   = pY0,                 # True counterfactual DRS 0
               pY1   = pY1,                 # True counterfactual DRS 1
               pY2   = pY2,                 # True counterfactual DRS 2
               Y     = Y)                   # Binary outcome
}

## Generate data given list of coefficients and N
GenerateData <- function(lstParams) {

    ## Generate a covariate matrix
    X          <- CreateX(lstParams$N)
    ## Assign treatment model output
    dfTrModel  <- AssignT(alpha10 = lstParams$alpha10, alpha20 = lstParams$alpha20,
                          alpha1X = lstParams$alpha1X, alpha2X = lstParams$alpha2X,
                          X = X)
    ## Assign outcome model output
    dfOutModel <- AssignY(beta0 = lstParams$beta0, betaX = lstParams$betaX,
                          betaT = lstParams$betaT, betaTX = lstParams$betaTX,
                          eMod = lstParams$eMod,
                          X = X, Tr = dfTrModel$Tr)
    ## Create a data frame containing covariates, treatment, and outcome
    df         <- as.data.frame(X)
    df         <- cbind(df, dfTrModel)
    df         <- cbind(df, dfOutModel)
    ## Return generated data with true PS and DRS
    df
}


## Scenario generator
## Give a list of named lists
## Each named list is a list of possible values (scalar or vector)
GenerateScenarios <- function(lstLstPossibleValues) {
    ## Scenario X Parameter data frame
    ## Each row represents a sceario
    dfScenarios <- do.call(expand.grid, lstLstPossibleValues)

    ## Number of scenarios
    nScenarios <- nrow(dfScenarios)

    ## Position of alphas (need special handling)
    posAlphas <- which(names(dfScenarios) == "alphas")

    ## Work on each scenario: Take i-th row and repackage into a list
    lstScenarios <- lapply(seq_len(nScenarios), function(i) {
        ## Each cell is a list, so these need c()ing
        ## This works for cells except the nested alphas
        lstExceptAlphas <- do.call(c, dfScenarios[i, -1 * posAlphas])
        ## alphas is nested, so it needs opening up before c()ing
        c(dfScenarios[i, posAlphas][[1]],
          lstExceptAlphas)
    })
}


## Function to generate R iterations and save for one scenario
GenerateRIterSave <- function(R, lstParams, scenarioCount) {

    ## Generate R iterations of the same scenario
    ## This is intentionally not parallelized for reproducibility
    lstData <- lapply(seq_len(R), function(i) {
        GenerateData(lstParams)
    })
    ## Generate file name
    fileName <- sprintf("Scenario%03d_R%d.RData",
                        scenarioCount, R)
    ## Save
    save(lstParams, scenarioCount, R, lstData, file = fileName)
    ## No return value
    NULL
}
