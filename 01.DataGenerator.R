#!/usr/local/bin/Rscript

################################################################################
### Generate three group data based on treatment model
##
## Created on: 2015-12-18
## Author: Kazuki Yoshida
################################################################################


### Prepare environment
################################################################################

## sink() if being run non-interactively
if (sink.number() != 0) {sink()}
..scriptFileName.. <- gsub("^--file=", "", Filter(function(x) {grepl("^--file=", x)}, commandArgs()))
if (length(..scriptFileName..) == 1) {
    sink(file = paste0(..scriptFileName.., ".txt"), split = TRUE)
    options(width = 100)
}

## Record start time
start_time <- Sys.time()
cat("### Started ", as.character(start_time), "\n")

## Configure parallelization
## Parallel backend for foreach (also loads foreach and parallel; includes doMC)
library(doParallel)
## Reproducible parallelization
library(doRNG)
## Detect core count
nCores <- min(parallel::detectCores(), 8)
## Used by parallel::mclapply() as default
options(mc.cores = nCores)
## Used by doParallel as default
options(cores = nCores)
## Register doParallel as the parallel backend for foreach
## http://stackoverflow.com/questions/28989855/the-difference-between-domc-and-doparallel-in-r
doParallel::registerDoParallel(cores = nCores)
## Report multicore use
cat("### Using", foreach::getDoParWorkers(), "cores\n")
cat("### Using", foreach::getDoParName(), "as backend\n")

## Load packages
library(magrittr)
library(dplyr)
library(tidyr)
library(tableone)


###
### Load data generating functions
################################################################################

source("./function_definitions/01.DataGeneratorFunctions.R")


###
### Scenarios for use
################################################################################

alphakX_base   <- c(1.0,1.0,0.1,2.0,2.0,0.4,0,0,0,0)
alphakX_nonlin <- c(0.8,0.06,0.06,1.5,1.5,0.4,0.8,0.04,0.08,1.0)
betaX_base   <- 1/2 * alphakX_base
betaX_nonlin <- 1/2 * alphakX_nonlin

### 33:33:33
alpha_cov1_E333333 <- list(XT_assoc = 1,
                           pT = c(`0` = 33, `1` = 33, `2` = 33),
                           alpha10 = -0.13,
                           alpha20 = -0.26,
                           alpha1X = 1/16 * alphakX_nonlin,
                           alpha2X = 1/8 * alphakX_nonlin)

alpha_cov2_E333333 <- list(XT_assoc = 2,
                           pT = c(`0` = 33, `1` = 33, `2` = 33),
                           alpha10 = -0.2,
                           alpha20 = -0.5,
                           alpha1X = 1/8 * alphakX_nonlin,
                           alpha2X = 1/4 * alphakX_nonlin)

alpha_cov3_E333333 <- list(XT_assoc = 3,
                           pT = c(`0` = 33, `1` = 33, `2` = 33),
                           alpha10 = -0.35,
                           alpha20 = -1.0,
                           alpha1X = 1/4 * alphakX_nonlin,
                           alpha2X = 1/2 * alphakX_nonlin)

alpha_cov4_E333333 <- list(XT_assoc = 4,
                           pT = c(`0` = 33, `1` = 33, `2` = 33),
                           alpha10 = -0.52,
                           alpha20 = -1.95,
                           alpha1X = 1/2 * alphakX_nonlin,
                           alpha2X = 1   * alphakX_nonlin)

alpha_cov5_E333333 <- list(XT_assoc = 5,
                           pT = c(`0` = 33, `1` = 33, `2` = 33),
                           alpha10 = -0.75,
                           alpha20 = -3.75,
                           alpha1X = 1 * alphakX_nonlin,
                           alpha2X = 2 * alphakX_nonlin)

### 10:45:45
alpha_cov1_E104545 <- list(XT_assoc = 1,
                           pT = c(`0` = 10, `1` = 45, `2` = 45),
                           alpha10 = 1.30,
                           alpha20 = 1.18,
                           alpha1X = 1/16 * alphakX_nonlin,
                           alpha2X = 1/8 * alphakX_nonlin)

alpha_cov2_E104545 <- list(XT_assoc = 2,
                           pT = c(`0` = 10, `1` = 45, `2` = 45),
                           alpha10 = 1.31,
                           alpha20 = 1.04,
                           alpha1X = 1/8 * alphakX_nonlin,
                           alpha2X = 1/4 * alphakX_nonlin)

alpha_cov3_E104545 <- list(XT_assoc = 3,
                           pT = c(`0` = 10, `1` = 45, `2` = 45),
                           alpha10 = 1.25,
                           alpha20 = 0.7,
                           alpha1X = 1/4 * alphakX_nonlin,
                           alpha2X = 1/2 * alphakX_nonlin)

alpha_cov4_E104545 <- list(XT_assoc = 4,
                           pT = c(`0` = 10, `1` = 45, `2` = 45),
                           alpha10 = 1.30,
                           alpha20 = 0.20,
                           alpha1X = 1/2 * alphakX_nonlin,
                           alpha2X = 1/1 * alphakX_nonlin)

alpha_cov5_E104545 <- list(XT_assoc = 5,
                           pT = c(`0` = 10, `1` = 45, `2` = 45),
                           alpha10 = 1.55,
                           alpha20 = -0.65,
                           alpha1X = 1 * alphakX_nonlin,
                           alpha2X = 2 * alphakX_nonlin)

### 10:10:80
alpha_cov1_E101080 <- list(XT_assoc = 1,
                           pT = c(`0` = 10, `1` = 10, `2` = 80),
                           alpha10 = -0.1,
                           alpha20 = 1.87,
                           alpha1X = 1/16 * alphakX_nonlin,
                           alpha2X = 1/8 * alphakX_nonlin)

alpha_cov2_E101080 <- list(XT_assoc = 2,
                           pT = c(`0` = 10, `1` = 10, `2` = 80),
                           lpha10 = -0.12,
                           alpha20 = 1.70,
                           alpha1X = 1/8 * alphakX_nonlin,
                           alpha2X = 1/4 * alphakX_nonlin)

alpha_cov3_E101080 <- list(XT_assoc = 3,
                           pT = c(`0` = 10, `1` = 10, `2` = 80),
                           alpha10 = -0.2,
                           alpha20 = 1.47,
                           alpha1X = 1/4 * alphakX_nonlin,
                           alpha2X = 1/2 * alphakX_nonlin)

alpha_cov4_E101080 <- list(XT_assoc = 4,
                           pT = c(`0` = 10, `1` = 10, `2` = 80),
                           alpha10 = 0.05,
                           alpha20 = 1.4,
                           alpha1X = 1/2 * alphakX_nonlin,
                           alpha2X = 1 * alphakX_nonlin)

alpha_cov5_E101080 <- list(XT_assoc = 5,
                           pT = c(`0` = 10, `1` = 10, `2` = 80),
                           alpha10 = 0.6,
                           alpha20 = 1.7,
                           alpha1X = 1 * alphakX_nonlin,
                           alpha2X = 2 * alphakX_nonlin)


###
### Experiment with parameter settings
################################################################################

if (FALSE) {

    lstParams <- c(alpha_cov5_E333333,
                   list(beta0   = log(0.20),
                        betaX   = 1/2.5*betaX_nonlin,
                        betaT   = c(0,0),
                        betaTX  = c(0,0),
                        eMod    = 4,
                        N       = 10^5))
    df <- GenerateData(lstParams)
    df$Y <- as.factor(df$Y)
    prop.table(table(df$Tr))
    library(tableone)
    print(CreateTableOne(vars = c(paste0("X",1:10),"Y","expLp","pY"),
                         strata = "Tr", data = df, test = FALSE), smd = TRUE)

}


###
### Generate scenarios
################################################################################

### Usable parameter settings
## Treatment model:
## {33:33:33, 10:45:45, 10:10:80} x 2 levels of confounder-tx assoc.
##
## Outcome model:
## beta0: {log(0.05), log(0.20)}
## betaX multiplication factor: fix at 1/2.5 as this is "biological"
## betaT: {(0,0), (log(0.9), log(0.6))} for null and beneficial scenarios
## betaTX (interaction): {0, log(0.8)} for null and additional benefit
##
## Estimation process:
## X1-X6 only and full X1-X10

## These part is a list of lists, so the nesting is deeper.
lstAlphas <- list(alpha_cov1_E333333, alpha_cov5_E333333,
                  alpha_cov1_E104545, alpha_cov5_E104545,
                  alpha_cov1_E101080, alpha_cov5_E101080)

## Create a list of named lists of possible values
## It has to be a list even if there is only one value
lstLstPossibleValues <- list(alphas  = lstAlphas,
                             beta0   = list(log(0.05), log(0.20)),
                             betaX   = list(1 / 2.5 * betaX_nonlin),
                             betaT   = list(c(0,0), c(log(0.9),log(0.6))),
                             betaTX  = list(c(0,0), c(log(0.7),log(0.5))),
                             eMod    = list(4),
                             N       = list(6*10^3))

## Generate list of scenarios
## Each element is a parameter set for a scenario
lstScenarios <- GenerateScenarios(lstLstPossibleValues)


###
### Generate datasets
################################################################################

nScenarios <- length(lstScenarios)
cat("###", nScenarios, "scenarios\n")

R <- 1000
cat("###", R, "iterations each\n")


## Start time date
dateStr <- format(start_time, "%Y%m%d")

## Create and move to a data folder
dirName <- paste0("./data/")
dir.create(path = dirName)
setwd(dirName)

## Generate scenarios requested R times each
## Parallelized using doRNG for reproducibility
set.seed(20151228*10)

junk <- foreach(i = seq_along(lstScenarios)) %dorng% {
    ## Generate data R times and save (return value NULL)
    cat("### Working on scenario:", i, "\n")
    GenerateRIterSave(R = R,
                      lstParams = lstScenarios[[i]],
                      scenarioCount = i)
}

## Record execution time
cat("### Started at:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("### Ended at  :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")


cat("###
### Peek first dataset in each scenario
################################################################################\n")

## File name pattern
filePat <- paste0("_R", R, ".RData")

## List files
dataFiles <- Filter(f = function(elt) {grepl(filePat, elt)}, x = dir(dirName))

cat("### Load first dataset from each data file\n")
lstData1 <- foreach(file = dataFiles) %dopar% {
    load(paste0(dirName, file))
    lstData[[1]]
}

cat("### Create a TableOne object for each dataset\n")
lstTabs <- foreach(df = lstData1) %dopar% {
    df$Y <- as.factor(df$Y)
    CreateTableOne(vars = c(paste0("X",1:10), "Y", "expLp", "pY"),
                   strata = "Tr", data = df, test = FALSE)
}

cat("### Print (do not parallelize)\n")
junk <- lapply(seq_along(lstTabs), function(i) {
    cat("\n### Scenario", i, "dataset 1 ### \n")
    print(lstTabs[[i]], smd = TRUE)
    print(lstScenarios[[i]])
})


################################################################################
cat("
###
### Record package versions etc
################################################################################\n")
print(sessionInfo())
## Record execution time
end_time <- Sys.time()
cat("### Started  ", as.character(start_time), "\n")
cat("### Finished ", as.character(end_time), "\n")
print(end_time - start_time)
## Stop sinking to a file if active
if (sink.number() != 0) {sink()}
