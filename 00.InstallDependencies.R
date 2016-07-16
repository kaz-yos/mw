#!/usr/local/bin/Rscript

################################################################################
### Install package dependencies required for simulation
##
## Created on: 2016-07-12
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

## Define installer
Inst <- function (PKG) {
    if (!(PKG %in% rownames(installed.packages()))) {
        cat("### Installing", PKG, "\n")
        install.packages(pkgs = PKG, dependencies = TRUE)
    }
}


cat("
###
### Installing dependencies
################################################################################\n")

Inst("VGAM")
Inst("boot")
Inst("doParallel")
Inst("doRNG")
Inst("dplyr")
Inst("ggplot2")
Inst("grid")
Inst("gridExtra")
Inst("magrittr")
Inst("rJava")
Inst("reshape2")
Inst("sandwich")
Inst("survey")
Inst("tableone")
Inst("tidyr")


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
