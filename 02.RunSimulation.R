#!/usr/local/bin/Rscript

################################################################################
### Run a simulation scenario given as an argument
##
## Created on: 2015-03-17
## Author: Kazuki Yoshida
################################################################################


### Prepare environment
################################################################################

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


###
### Load necessary components
################################################################################
source("./function_definitions/07.Simulate.R")


###
### Invoke simulation with a data file name
################################################################################

## Data file name should be given as the first argument for this script
dataFileName <- commandArgs(trailingOnly = TRUE)[1]
## Stop if not data file name was given
stopifnot(!is.na(dataFileName))

## Load data
load(dataFileName)

## Iterate over the list of data frames
## Parallelized using mclapply()
dfOut <- Iterate(lstData = lstData,
                 psFormula = Tr ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                 lstParams = lstParams)


###
### Save all objects
################################################################################

save(dfOut, lstParams, R, scenarioCount,
     file = gsub(pattern = "Scenario", replacement = "Result", x = dataFileName))


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
