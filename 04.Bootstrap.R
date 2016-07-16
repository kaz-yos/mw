#!/usr/local/bin/Rscript

################################################################################
### Bootstrap for variance
##
## Created on: 2016-01-19
## Author: Kazuki Yoshida
################################################################################


### Prepare environment
################################################################################

startTime <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
cat("### Started at:", startTime, "\n")


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
## Register doParallel as the parallel backend with foreach
## http://stackoverflow.com/questions/28989855/the-difference-between-domc-and-doparallel-in-r
doParallel::registerDoParallel(cores = nCores)
## Report multicore use
cat("### Using", foreach::getDoParWorkers(), "cores\n")
cat("### Using", foreach::getDoParName(), "as backend\n")


## Load packages

## Load functions
source("./function_definitions/09.BootstrapFunctions.R")


## Configure sink()
if (sink.number() != 0) {sink()}
..scriptFileName.. <- gsub("^--file=", "", Filter(function(x) {grepl("^--file=", x)}, commandArgs()))
if (length(..scriptFileName..) == 1) {
    sink(file = paste0(..scriptFileName.., ".txt"), split = TRUE)
}
options(width = 120)



cat("###
### Load data
################################################################################\n")

## Argument 1:
dataFileName <- commandArgs(trailingOnly = TRUE)[1]
## Stop if not data file name was given
stopifnot(!is.na(dataFileName))

## Argument 2: Which part out of 1:10
part10 <- as.integer(commandArgs(trailingOnly = TRUE)[2])
## Must be NA or 1:10
stopifnot(is.na(part10) | part10 %in% c(1:10))


## Load data
load(dataFileName)

cat("### Parameters\n")
lstParams


## Create index
if (!is.na(part10)) {
    cat("### Subset data: part ", part10, "of 10\n")
    ## Assumes R is in data file
    partSize <- R/10

    ## 1:partSize
    i_init <- seq_len(partSize)

    if (part10 == 10) {
        ## If it is the last part, include up to the very last dataset
        i_final <- seq(from = (partSize * 9) + 1, to = R, by = 1)
    } else {
        ## If not go by the partSize
        i_final <- i_init + partSize * (part10 - 1)
    }
    cat("### Working on datasets", min(i_final), "through", max(i_final), "\n")
} else {
    cat("### Working on all datasets\n")
    i_final <- seq_len(R)
}


cat("###
### Perform bootstrapping
################################################################################\n")

## Bootstrap iterations for each dataset
B <- 1000
cat("### ", B, "boostrap iterations\n")

## Create a closure
BootFun <- ConstructBootFun(psFormula = Tr ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                            lstParams = lstParams)

## Parallel execution
set.seed(201601190)
lstBootRes <- foreach(i = i_final) %dorng% {
    cat("### Running dataset", i, "\n")
    bootOut <- boot(data = lstData[[i]],
                    statistic = BootFun,
                    R = B, parallel = "no")
    ## B x p matrix of bootstrapped coefficients
    bootOut$t
}

cat("### Show first dataset results\n")
lstBootRes[[1]]


cat("###
### Save
################################################################################\n")

## Create save file name
if (!is.na(part10)) {
    repl <- sprintf("_part%02d.RData", part10)
    saveFileName <- gsub(pattern = "Scenario", replacement = "Bootstrap", x = dataFileName) %>%
        gsub(pattern = ".RData$", replacement = repl, x = .)
} else {
    saveFileName <- gsub(pattern = "Scenario", replacement = "Bootstrap", x = dataFileName)
}


dateString <- system("date '+%Y%m%d_%H%M%S'", intern = TRUE)
save(lstBootRes, lstParams, R, scenarioCount,
     file = saveFileName)


cat("### Started at:", startTime, "\n")
cat("### Ended at  :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

system("uptime")



################################################################################
cat("\n### Record package versions\n")
print(sessionInfo())
## Stop sinking to a file if active
if (sink.number() != 0) {sink()}
