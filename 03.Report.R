#!/usr/local/bin/Rscript

################################################################################
### Format files given as arguments and report results from simulation
##
## Created on: 2015-04-09
## Author: Kazuki Yoshida
################################################################################


### Prepare environment
################################################################################

## Configure sink()
if (sink.number() != 0) {sink()}
..scriptFileName.. <- gsub("^--file=", "", Filter(function(x) {grepl("^--file=", x)}, commandArgs()))
if (length(..scriptFileName..) == 1) {
    sink(file = paste0(..scriptFileName.., ".txt"), split = TRUE)
    options(width = 120)
}

## Record start time
start_time <- Sys.time()
cat("### Started ", as.character(start_time), "\n")

## Load packages
library(magrittr)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)

source("./function_definitions/05.BalanceCheck.R")
source("./function_definitions/07.Simulate.R")
source("./function_definitions/08.ReportingFunctions.R")


cat("###
### Load data files
################################################################################\n")

## List all .data/Result* files
resFiles <- system("ls ./data/Result*.RData", intern = TRUE)


## Stop if no file is specified
stopifnot(length(resFiles) > 0 & all(!is.na(resFiles)))

cat("### Files to be read\n")
print(resFiles)

## Load into a list
lstRes <- lapply(resFiles, function(file) {
    load(file)
    list(dfOut = dfOut, lstParams = lstParams, R = R, scenarioCount = scenarioCount)
})

reports <- lapply(lstRes, function(res) {
    Report(scenarioCount = res$scenarioCount,
           lstParams     = res$lstParams,
           dfOut         = res$dfOut)
}) %>% do.call(rbind, .) %>% as.data.frame


cat("###
### Name scenarios
################################################################################\n")

## Original sample size
reports$N <- factor(reports$N)

## Strength of covariate-treatment association (stronger = poor covariate overlap)
reports$XT_assoc <- factor(reports$XT_assoc, levels = c(1,5),
                           labels = c("Good overlap", "Poor overlap"))

## Exposure prevalence
reports$pExpo <- factor(reports$pT.2, levels = c(33,45,80),
                        labels = c("33:33:33", "10:45:45", "10:10:80"))

## Main effects
reports$effects <- factor((reports$betaT1 != 0), levels = c(FALSE,TRUE),
                          labels = c("Null main effects", "Non-null main effects"))

## Effect modification
reports$em <- factor((reports$betaTX1 != 0), levels = c(FALSE,TRUE),
                     labels = c("Modification (-)", "Modification (+)"))

## Disease prevalence
reports$pDis <- factor(reports$beta0, levels = sort(unique(reports$beta0)),
                       labels = exp(sort(unique(reports$beta0))))
## Combine with N
## reports$pDis <- interaction(reports$pDis, reports$N, sep = ":")

## Specify index variables
indexVars <- c("scenario","N","XT_assoc","pExpo","effects","em","pDis")

cat("### Show scenarios\n")
reports[indexVars]
save(reports, indexVars, file = "./data/Report.RData")


cat("###
### Numerical Examination
################################################################################\n")

cat("### Check if there are missing values (corrupt scenarios)\n")
colSumNA <- colSums(is.na(reports))
colSumNA[colSumNA != 0]


cat("### Check if caliper widening ever happened\n")
reports[reports$caliperPw > 0,c("scenario","XT_assoc",paste0("pT.", c(0,1,2)))]


cat("### Magnitude of bias\n")
## Better if closer to 1 (closer to 0 on log scale)
absLogBiasM <- abs(log(reports[,PickMatchNames(names(reports), "M.biasRR")]))
absLogBiasMw <- abs(log(reports[,PickMatchNames(names(reports), "Mw.biasRR")]))
absLogBiasRatio <- absLogBiasM / absLogBiasMw
cat("### Proportion M > Mw using entire sample true value\n")
mean(absLogBiasRatio[,1:3] > 1)
cat("### Proportion M > Mw using common support true value\n")
mean(absLogBiasRatio[,4:6] > 1)


cat("### Magnitude of variance\n")
## Better if closer to 1 (closer to 0 on log scale)
trueVarM  <- reports[,PickMatchNames(names(reports), "M.trueV")]
trueVarMw <- reports[,PickMatchNames(names(reports), "Mw.trueV")]
trueVarIp <- reports[,PickMatchNames(names(reports), "Ip.trueV")]
cat("### Proportion M > Mw\n")
mean((trueVarM / trueVarMw)[,1:3] > 1)
cat("### Proportion Ip > Mw\n")
mean((trueVarIp / trueVarMw)[,1:3] > 1)
cat("### Proportion M > Ip\n")
mean((trueVarM / trueVarIp)[,1:3] > 1)


cat("### Magnitude of MSE\n")
## Better if closer to 1 (closer to 0 on log scale)
mseM  <- reports[,PickMatchNames(names(reports), "M.mse\\.")]
mseMw <- reports[,PickMatchNames(names(reports), "Mw.mse\\.")]
mseIp <- reports[,PickMatchNames(names(reports), "Ip.mse\\.")]
cat("### Proportion M > Mw\n")
mean((mseM / mseMw)[,1:3] > 1)
cat("### Proportion Ip > Mw\n")
mean((mseIp / mseMw)[,1:3] > 1)
cat("### Proportion M > Ip\n")
mean((mseM / mseIp)[,1:3] > 1)


cat("### False positive rate\n")
## Better if closer to 1 (closer to 0 on log scale)
dat <- Gather(reports, indexVars, "pRej")
dat <- subset(dat, effects == "Null main effects" & em == "Modification (-)")
dat %>%
    group_by(method) %>%
    summarize(mean      = mean(value),
              totalNull = n(),
              antiCon05 = sum(value > 0.05),
              antiCon06 = sum(value > 0.06),
              antiCon07 = sum(value > 0.07))



cat("###
### Graphical Examination
################################################################################\n")

## Name graph using the date and time
pdf(file = "./figures/_ExperimentalFigures.pdf", width = 14, height = 10, family = "sans")

## Graph prototype
gg <- ggplot(mapping = aes(x = key, y = value)) +
    geom_point() +
    labs(y = "", x = "") +
    theme_bw() + theme(legend.key = element_blank(),
                       axis.text.x = element_text(angle = 90))


cat("### Sample size\n")
## dat <- Gather(reports, indexVars, "\\.n$")
## gg %+% dat +
##     aes(group = scenario, color = pExpo, x = method) +
##     geom_line(size = 0.1) +
##     facet_grid(. ~ XT_assoc) +
##     labs(title = "Sample size") +
##     geom_hline(yintercept = reports$U.n[1])

## dat <- Gather(reports, indexVars, "\\.nCs$")
## gg %+% dat +
##     aes(group = scenario, color = pExpo, x = method) +
##     geom_line(size = 0.1) +
##     facet_grid(. ~ XT_assoc) +
##     labs(title = "Sample size (common support)") +
##     geom_hline(yintercept = reports$U.n[1])

dat <- Gather(reports, indexVars, "\\.n$|\\.nCs$")
dat$key <- factor(gsub(".*\\.", "", as.character(dat$key)), levels = c("n","nCs"),
                  labels = c("All", "Within common support"))
gg %+% dat +
    aes(group = scenario, color = pExpo, x = method) +
    geom_line(size = 0.1) +
    facet_grid(XT_assoc ~ key) +
    labs(title = "Sample size (All vs Within common support)") +
    scale_y_continuous(limit = c(0,NA)) +
    geom_hline(yintercept = unique(reports$U.n))


cat("### SMD\n")
dat <- Gather(reports, indexVars, "smd")
dat$key <- gsub("^.*\\.", "", as.character(dat$key))
dat$key <- factor(dat$key, levels = paste0("X", 1:10))
gg %+% dat +
    labs(title = "SMD") +
    aes(group = scenario, color = pExpo) +
    geom_line() +
    facet_grid(XT_assoc ~ method) +
    scale_y_continuous() +
    coord_cartesian(ylim = c(0, 0.5)) +
    geom_hline(yintercept = 0.10, size = 0.01) +
    geom_hline(yintercept = 0.00)


cat("### Prevalence of modifier\n")
dat <- Gather(reports, indexVars, "pMod$")
gg %+% dat +
    aes(group = scenario, color = pExpo, x = method) +
    geom_line(alpha = 0.5) +
    scale_y_continuous(limit = c(0, NA)) +
    facet_grid(. ~ XT_assoc) +
    labs(title = "Prevalence of modifier", x = "")



cat("### Bias related\n")
dat <- Gather(reports, indexVars, "coef")
dat$value <- exp(dat$value)
dat$trueValue <- Gather(reports, indexVars, "trueRR\\.")$value
gg %+% dat +
    aes(group = scenario, color = pExpo, shape = pDis, x = method) +
    labs(x = "") +
    ## To indicate the true values
    geom_point(mapping = aes(y = trueValue), shape = 4, alpha = 1/5, size = 5) +
    geom_line() +
    scale_y_log10(breaks = c(1/4, 1/2, 1, 2, 4)) +
    coord_cartesian(ylim = c(1/4, 4)) +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    labs(title = "RR along with with true RR", x = "")

dat <- Gather(reports, indexVars, "trueRR\\.")
gg %+% dat +
    aes(group = scenario, color = pExpo, shape = pDis, x = method) +
    geom_line() +
    scale_y_log10(breaks = c(1/4, 1/2, 1, 2, 4)) +
    coord_cartesian(ylim = c(1/4, 4)) +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    labs(title = "True RR", x = "")

dat <- Gather(reports, indexVars, "biasRR\\.")
gg %+% dat +
    aes(group = scenario, color = pExpo, shape = pDis, x = method) +
    geom_line() +
    scale_y_log10(breaks = c(1/4, 1/2, 1, 2, 4)) +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    coord_cartesian(ylim = c(1/4, 4)) +
    geom_hline(yintercept = 1) +
    labs(title = "Bias (RR / true RR)")



cat("### Variance\n")
dat <- Gather(reports, indexVars, "trueV")
gg %+% dat +
    aes(group = scenario, color = pExpo, shape = pDis, x = method) +
    geom_line() +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    coord_cartesian(ylim = c(0, 2)) +
    labs(title = "True variance", x = "")

dat <- Gather(reports, indexVars, "vars")
gg %+% dat +
    aes(group = scenario, color = pExpo, shape = pDis, x = method) +
    geom_line() +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    coord_cartesian(ylim = c(0, 2)) +
    labs(title = "Mean estimated variance", x = "")

datWide <-
    reports[, c(PickMatchNames(names(reports), "vars"))] /
    reports[, c(PickMatchNames(names(reports), "trueV"))]
datWide <- cbind(datWide, reports[, indexVars])
dat <- Gather(datWide, indexVars, "vars")
gg %+% dat +
    aes(group = scenario, color = pExpo, shape = pDis, x = method) +
    geom_line() +
    facet_grid(effects + XT_assoc ~ em + contrast) +
        geom_hline(yintercept = 1) +
    labs(title = "Mean estimated variance / true variance", x = "")


cat("### MSE\n")
dat <- Gather(reports, indexVars, "mse\\.")
gg %+% dat +
    aes(group = scenario, color = pExpo, shape = pDis, x = method) +
    geom_line() +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    coord_cartesian(ylim = c(0, 1.5)) +
    labs(title = "MSE", x = "")


cat("### Alpha rate\n")
dat <- Gather(reports, indexVars, "pRej")
dat <- subset(dat, effects == "Null main effects" & em == "Modification (-)")
gg %+% dat +
    aes(group = scenario, color = pExpo, shape = pDis, x = method) +
    geom_line() +
    facet_grid(XT_assoc ~ contrast) +
    scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
    ## coord_cartesian(ylim = c(0, 0.25)) +
    geom_hline(yintercept = 0.05) +
    labs(title = "Observed alpha error rate (null scenarios)", x = "")


cat("### Coverage\n")
dat <- Gather(reports, indexVars, "cvr\\.")
gg %+% dat +
    aes(group = scenario, color = pExpo, shape = pDis, x = method) +
    geom_line() +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    ## coord_cartesian(ylim = c(0.75, 1.00)) +
    geom_hline(yintercept = 0.95) +
    labs(title = "Coverage", x = "")


cat("### Minimum counts\n")
dat <- Gather(reports, indexVars, "NCases\\.")
dat$key <- gsub("minNCases\\.", "", as.character(dat$key))
gg %+% dat +
    aes(group = scenario, color = pExpo, shape = pDis, x = key) +
    geom_line() +
    facet_grid(effects + XT_assoc ~ em + pDis) +
    coord_cartesian(ylim = c(0, 25)) +
    geom_hline(yintercept = 0) +
    labs(title = "Minimum case counts", x = "")

dat <- Gather(reports, indexVars, "minX10\\.")
dat$key <- gsub("minX10\\.", "", as.character(dat$key))
gg %+% dat +
    aes(group = scenario, color = pExpo, shape = pDis, x = key) +
    geom_line() +
    facet_grid(effects + XT_assoc ~ em + pDis) +
    coord_cartesian(ylim = c(0, 25)) +
    geom_hline(yintercept = 0) +
    labs(title = "Minimum X10 counts", x = "")

dev.off()


cat("
###
### Production figures
################################################################################\n")

## http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
## Color-blind friendly palette with black.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(4,6,7)]


gg <- ggplot(mapping = aes(x = key, y = value)) +
    labs(title = NULL, y = NULL, x = NULL) +
    scale_color_manual(values = cbbPalette) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme_bw() + theme(legend.key = element_blank(),
                       legend.key.width = unit(1, "cm"),
                       legend.position = "bottom",
                       strip.background = element_blank(),
                       panel.margin = unit(0.5, "lines"))

###  Figures

###   Average sample size
pdf(file = "./figures/Figure2.pdf", width = 5, height = 5, family = "sans")
dat <- Gather(reports, indexVars, "\\.n$")
## Average across random variability with essentially same scenarios
dat <- dat %>%
    group_by(pExpo,method,key,XT_assoc) %>%
    summarize(value = mean(value))
gg %+% dat +
    aes(group = pExpo, linetype = pExpo, x = method) +
    geom_line() +
    geom_point() +
    facet_grid(. ~ XT_assoc) +
    geom_hline(yintercept = c(reports$U.n[1], 0), alpha = 1/5) +
    theme(panel.margin = unit(0, "lines"),
          legend.margin = unit(0, "lines")) +
    labs(title = NULL, x = NULL, y = "Sample Size")
dev.off()


###   SMD
pdf(file = "./figures/Figure3.pdf", width = 5, height = 5, family = "sans")
dat <- Gather(reports, indexVars, "smd")
dat$key <- gsub("^.*\\.", "", as.character(dat$key))
dat$key <- factor(dat$key, levels = paste0("X", 1:10))
## Restrict to 3 covariates
dat <- subset(dat, key %in% paste0("X", c(1,4,7)))
## Average across random variability with essentially same scenarios
dat <- dat %>%
    group_by(pExpo,method,key,XT_assoc) %>%
    summarize(value = mean(value))
gg %+% dat +
    aes(group = pExpo, linetype = pExpo, x = method) +
    geom_line(alpha = 2/3) +
    geom_point() +
    facet_grid(XT_assoc ~ key, scales = "free") +
    scale_y_continuous() +
    coord_cartesian(ylim = c(0, 0.5)) +
    geom_hline(yintercept = 0.10, size = 0.3, alpha = 3/5) +
    geom_hline(yintercept = 0.00, alpha = 3/5) +
    labs(title = NULL, x = NULL, y = "Absolute Standardized Mean Difference")
dev.off()


###  eFigures

###   Bias
pdf(file = "./figures/eFigure2.pdf", width = 10, height = 10, family = "sans")
dat <- Gather(reports, indexVars, "biasRR\\.")
gg %+% dat +
    aes(group = scenario, linetype = pExpo, shape = pDis, x = method) +
    geom_line(alpha = 2/3) +
    geom_point(alpha = 2/3) +
    scale_y_log10(breaks = c(1/4, 1/2, 3/4, 1, 1.5, 2, 3, 4)) +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    coord_cartesian(ylim = c(3/4, 3)) +
    geom_hline(yintercept = 1, alpha = 1/5) +
    labs(title = NULL, x = NULL, y = "Bias (Estimated Risk Ratio / True Risk Ratio)")
dev.off()


###   True RR
pdf(file = "./figures/eFigure3.pdf", width = 10, height = 10, family = "sans")
dat <- Gather(reports, indexVars, "trueRR\\.")
gg %+% dat +
    aes(group = scenario, linetype = pExpo, x = method) +
    geom_line(alpha = 2/3) +
    scale_y_log10(breaks = c(1/4, 2/5, 1/2, 3/4,  1, 2, 4)) +
    coord_cartesian(ylim = c(exp(-1), 1.1)) +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    labs(title = NULL, x = NULL, y = "True Risk Ratio")
dev.off()


###   True variance
pdf(file = "./figures/eFigure4.pdf", width = 10, height = 10, family = "sans")
dat <- Gather(reports, indexVars, "trueV")
gg %+% dat +
    aes(group = scenario, linetype = pExpo, shape = pDis, x = method) +
    geom_line() +
    geom_point() +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = NULL, x = NULL, y = "True Variance")
dev.off()


###   Mean estimated variance
pdf(file = "./figures/eFigure5.pdf", width = 10, height = 10, family = "sans")
dat <- Gather(reports, indexVars, "vars")
gg %+% dat +
    aes(group = scenario, linetype = pExpo, shape = pDis, x = method) +
    geom_line() +
    geom_point() +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = NULL, x = NULL, y = "Estimated Variance")
dev.off()

###   eFigure6.pdf is bootstrap variance figure created elsewhere


###   MSE
pdf(file = "./figures/eFigure7.pdf", width = 10, height = 10, family = "sans")
dat <- Gather(reports, indexVars, "mse\\.")
gg %+% dat +
    aes(group = scenario, linetype = pExpo, shape = pDis, x = method) +
    geom_line() +
    geom_point() +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    coord_cartesian(ylim = c(0, 1.5)) +
    labs(title = NULL, x = NULL, y = "Mean Squared Error")
dev.off()


###   Type I error rate
pdf(file = "./figures/eFigure8.pdf", width = 5, height = 5, family = "sans")
dat <- Gather(reports, indexVars, "pRej")
dat <- subset(dat, effects == "Null main effects" & em == "Modification (-)")
gg %+% dat +
    aes(group = scenario, linetype = pExpo, shape = pDis, x = method) +
    geom_line() +
    geom_point() +
    facet_grid(XT_assoc ~ contrast) +
    scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
    ## coord_cartesian(ylim = c(0, 0.25)) +
    geom_hline(yintercept = 0.05) +
    labs(title = NULL, x = NULL, y = "Type I Error Rate")
dev.off()


###   Coverage probability
pdf(file = "./figures/eFigure9.pdf", width = 10, height = 10, family = "sans")
dat <- Gather(reports, indexVars, "cvr\\.")
gg %+% dat +
    aes(group = scenario, linetype = pExpo, shape = pDis, x = method) +
    geom_line() +
    geom_point() +
    facet_grid(effects + XT_assoc ~ em + contrast) +
    ## coord_cartesian(ylim = c(0.75, 1.00)) +
    geom_hline(yintercept = 0.95) +
    labs(title = NULL, x = NULL, y = "Coverage Probability")
dev.off()

###   eFigure10.pdf is empirical SMD figure created elsewhere


cat("
###
### Anomaly assessment
################################################################################\n")

### High variance for M
cat("### Assessment of very high variance in matching\n")
cat("### Show configuration for most extreme scenario\n")
scenarioNum <- which.max(reports$M.trueV.2v0)
reports[scenarioNum, c("pExpo","XT_assoc","effects","em","pDis")]

cat("### Extract iteration results\n")
res <- lstRes[[scenarioNum]]

cat("### Show coefficients distribution (very low coefficients)\n")
summary(res$dfOut[,PickMatchNames(names(reports), "coef")])

cat("### Show most extreme iterations\n")
res$dfOut[res$dfOut$M.coef.2v0 < min(res$dfOut$M.coef.2v0) + 1.0, ]

cat("### Iteration count of most extreme iteration\n")
as.numeric(rownames(res$dfOut[which.min(res$dfOut$M.coef.2v0),]))


### NA's in performance assessment (SMD for X10)
cat("### Check columns having NA's (This should be empty)\n")
countNa <- foreach(lst = lstRes, .combine = rbind) %do% {
        colSums(is.na(lst$dfOut))
}
countNa[, colSums(countNa > 0) > 0]


### False positive scenarios
cat("### Check which scenarios are producing false positives\n")

reportsNull <- subset(reports, effects == "Null main effects" & em == "Modification (-)")

FalsePosScenarios <- function(indexVars, data, pat) {
    vars <- PickMatchNames(names(data), pat)
    list(alpha005 = data[rowSums(data[,vars] > 0.05) > 0, c(indexVars, vars)],
         alpha006 = data[rowSums(data[,vars] > 0.06) > 0, c(indexVars, vars)],
         alpha007 = data[rowSums(data[,vars] > 0.07) > 0, c(indexVars, vars)])
}

cat("### Violating scenarios for M\n")
FalsePosScenarios(indexVars, reportsNull, "M.pRej")

cat("### Violating scenarios for Mw\n")
FalsePosScenarios(indexVars, reportsNull, "Mw.pRej")


### Undercoverage scenarios
cat("### Check which scenarios are producing undercoverage\n")
UnderCiScenarios <- function(indexVars, data, pat) {
    vars <- PickMatchNames(names(data), pat)
    ## All over coverage
    list(cvr_all_gr097 = data[rowSums(data[,vars] > 0.96) > 2, c(indexVars, vars)],
         cvr_all_gr096 = data[rowSums(data[,vars] > 0.96) > 2, c(indexVars, vars)],
         ## Any undercoverage
         cvr_ls095 = data[rowSums(data[,vars] < 0.95) > 0, c(indexVars, vars)],
         cvr_ls094 = data[rowSums(data[,vars] < 0.94) > 0, c(indexVars, vars)],
         cvr_ls093 = data[rowSums(data[,vars] < 0.93) > 0, c(indexVars, vars)])
}

cat("### Violating scenarios for M\n")
UnderCiScenarios(indexVars, reports, "M.cvr\\.")

cat("### Violating scenarios for Mw\n")
UnderCiScenarios(indexVars, reports, "Mw.cvr\\.")


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
