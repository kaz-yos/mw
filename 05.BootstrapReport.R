#!/usr/local/bin/Rscript

################################################################################
### Script to report bootstrap results
##
## Created on: 2016-01-21
## Author: Kazuki Yoshida
################################################################################


### Prepare environment
################################################################################

## Configure sink()
if (sink.number() != 0) {sink()}
..scriptFileName.. <- gsub("^--file=", "", Filter(function(x) {grepl("^--file=", x)}, commandArgs()))
if (length(..scriptFileName..) == 1) {
    sink(file = paste0(..scriptFileName.., ".txt"), split = TRUE)
}
options(width = 120)

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
## Register doParallel as the parallel backend with foreach
## http://stackoverflow.com/questions/28989855/the-difference-between-domc-and-doparallel-in-r
doParallel::registerDoParallel(cores = nCores)
## Report multicore use
cat("### Using", foreach::getDoParWorkers(), "cores\n")
cat("### Using", foreach::getDoParName(), "as backend\n")


## Load packages
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

## Load functions
source("./function_definitions/08.ReportingFunctions.R")


cat("
###
### Load data
################################################################################\n")

## Load all ./data/Bootstrap* files
fileNames <- system("ls ./data/Bootstrap*.RData", intern = TRUE)

## Keep part-by-part bootstrap files only
bootFileNamesDf <- data.frame(fileName = Filter(f = function(elt) {grepl("_part", elt)}, x = fileNames), stringsAsFactors = FALSE)

## Parsers
ScenarioNumber <- function(fileName) {
    fileName %>%
        gsub(".*Bootstrap", "", .) %>%
        gsub("_.*", "", .) %>%
        as.numeric
}
SegmentNumber <- function(fileName) {
    fileName %>%
        gsub("^.*_part", "", .) %>%
        gsub("\\.RData", "", .) %>%
        as.numeric
}
bootFileNamesDf$scenario <- ScenarioNumber(bootFileNamesDf$fileName)
bootFileNamesDf$segment  <- SegmentNumber(bootFileNamesDf$fileName)


### Load existing files
N <- 48

lstLstBootRes <- lapply(bootFileNamesDf$fileName, function(fileName) {
    load(fileName)
    ## list of 100 matrices (each has 1000 bootstrap iterations * 3 parameters)
    lstBootRes
})

## Create mean bootstrap variance at each 1/10-th segment
dfBootVars <- lapply(lstLstBootRes, function(lst) {
    lapply(lst, function(mat) {
        apply(mat, MARGIN = 2, FUN = var)
    }) %>%
        do.call(rbind, .) %>%
        colMeans(.)
}) %>%
    do.call(rbind, .) %>%
    as.data.frame
names(dfBootVars) <- paste0("Mw.bootV.", c("1v0","2v0","2v1"))

## Create mean bootstrap mean at each 1/10-th segment
dfBootMeans <- lapply(lstLstBootRes, function(lst) {
    lapply(lst, function(mat) {
        apply(mat, MARGIN = 2, FUN = mean)
    }) %>%
        do.call(rbind, .) %>%
        colMeans(.)
}) %>%
    do.call(rbind, .) %>%
    as.data.frame
names(dfBootMeans) <- paste0("Mw.coefB.", c("1v0","2v0","2v1"))


## Finalize report dataset
reportsBoot <- cbind(bootFileNamesDf[,-1], dfBootVars, dfBootMeans) %>%
    plyr::ddply(.data = .,
                .variables = "scenario",
                .fun = function(df) {
                    colMeans(df[c(names(dfBootVars), names(dfBootMeans))])
                })




### Load existing data
load("./data/Report.RData")

### inner join
reportsMerge <- full_join(reports, reportsBoot, by = "scenario")


## Gather columns for variance
reportsVars <- Gather(reportsMerge, indexVars, "Mw\\.vars|Mw.*V")
## Transform keys
reportsVars$key <- as.character(reportsVars$key) %>%
    gsub("Mw.", "", .) %>%
    gsub("\\..*", "", .)
## Give appropriate ordering
reportsVars$key <-
    factor(reportsVars$key, levels = c("vars", "trueV", "bootV"),
           labels = c("Est.", "True", "Boot."))


## Gather columns for mean
reportsCoefs <- Gather(reportsMerge, indexVars, "Mw\\.trueCoef\\.|Mw\\.coef\\.|Mw\\.coefB")
## Transform keys
reportsCoefs$key <- as.character(reportsCoefs$key) %>%
    gsub("Mw.", "", .) %>%
    gsub("\\..*", "", .)
## Give appropriate ordering
reportsCoefs$key <-
    factor(reportsCoefs$key, levels = c("coef", "trueCoef", "coefB"),
           labels = c("Est.", "True", "Boot."))


cat("###
### Summarize graphically
################################################################################\n")

## http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
## Color-blind friendly palette with black.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(4,6,7)]

gg <- ggplot(mapping = aes(x = key, y = value)) +
    labs(y = NULL, x = NULL) +
    scale_color_manual(values = cbbPalette) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme_bw() + theme(legend.key = element_blank(),
                       legend.key.width = unit(1, "cm"),
                       legend.position = "bottom",
                       strip.background = element_blank(),
                       panel.margin = unit(0.5, "lines"))


## Variance figure
pdf(file = "./figures/eFigure5.pdf", width = 10, height = 10, family = "sans")
gg %+% reportsVars +
    aes(group = scenario, linetype = pExpo, shape = pDis, x = key) +
    geom_line() +
    geom_point() +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    ## coord_cartesian(ylim = c(0, 2)) +
    labs(title = NULL, y = "Variance", x = NULL)
dev.off()

## Corresponding average coefficient figure for sanity check
pdf(file = "./figures/eFigure5_SanityCheck.pdf", width = 10, height = 10, family = "sans")
gg %+% reportsCoefs +
    aes(group = scenario, linetype = pExpo, shape = pDis, x = key) +
    geom_line() +
    geom_point() +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    ## coord_cartesian(ylim = c(0, 2)) +
    labs(title = "Bootstrap Mean Comparison", x = "")
dev.off()


################################################################################
cat("\n### Record package versions\n")
print(sessionInfo())
## Record execution time
end_time <- Sys.time()
cat("### Started  ", as.character(start_time), "\n")
cat("### Finished ", as.character(end_time), "\n")
print(end_time - start_time)
## Stop sinking to a file if active
if (sink.number() != 0) {sink()}
