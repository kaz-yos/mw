################################################################################
### R Frontend Function for Rassen 2013 Matching Algorithm
## Just an interface between R and Rassen 2013 Java backend jar file
## Created on: 2015-03-16
## Author: Kazuki Yoshida
################################################################################


### Necessary packages

library(magrittr)


### Problem in parallelization
## problem using rJava with parallel::mclapply
## http://r.789695.n4.nabble.com/problem-using-rJava-with-parallel-mclapply-td4680245.html
## Running mclapply with rJava: Load rJava within each fork
## https://hpctech.wordpress.com/2014/05/28/running-mclapply-with-rjava/

Match3WayRassen <- function(cpath = "./rassen_toolbox/java/pharmacoepi.jar",
                            matchType = "nn", cohort,
                            numGroups = 3, groupNames = c(0,1,2),
                            matchRatio = 3, caliper = 0.05,
                            startDigit = 5, endDigit = 1) {

    ## rJava should be loaded here within each fork of parallelization
    library(rJava)

    ## Initialize Java VM
    ## What does Java option -Xmx stand for?
    ## http://stackoverflow.com/questions/5374455/what-does-java-option-xmx-stand-for
    .jinit(classpath = cpath, parameters = c("-Xmx8000m"), force.init = FALSE)

    ## Java null object reference
    m <- .jnull("org.drugepi.match.Match")

    ## Create a Java object (Instantiate match object)
    m <- .jnew("org.drugepi.match.Match")

    ## print(sprintf("Match version %s", m$version))

    ## Call a Java method, initMatch
    ## https://github.com/bwh-dope/pharmacoepi_toolbox/blob/eecf267b3c95410a6b6aedb8a55c8540d7ea7e43/src/org/drugepi/match/Match.java#L127
    .jcall(obj = m, returnSig = "V", method = "initMatch",
           ## parameters passed to the Java method
           ## String matchType, String numGroups
           matchType, as.character(numGroups))

    ## Set configs
    ## Obtains the value of a field
    .jfield(m, "matchRatio")           <- as.integer(matchRatio)
    .jfield(m, "fixedRatio")           <- as.integer(1)
    .jfield(m, "parallelMatchingMode") <- as.integer(1)
    ## not all options will be used in all runs, but set them anyway
    .jfield(m, "caliper")              <- as.double(caliper)
    .jfield(m, "startDigit")           <- as.integer(startDigit)
    .jfield(m, "endDigit")             <- as.integer(endDigit)

    ## treatment groups
    for (groupName in groupNames) {
        ## Last one will be reference
        .jcall(m, "V", "addMatchGroup", as.character(groupName))
    }

    ## Set output file
    ## Add additional random number to avoid hitting upper limit of names
    pattern <- paste0("pe_toolbox_", round(runif(n = 1) * 10^6), "_")
    outfilePath <- tempfile(pattern = pattern)
    .jfield(m, "outfilePath") <- outfilePath

    ## Convert R data frame to a huge string
    ## ID, Treatment group, PS1, (PS2)
    match_header <- paste(names(cohort), collapse = "\t")
    match_data   <- paste(do.call(function(...) {paste(sep = "\t", ...)}, cohort),
                          collapse="\n")

    ##
    .jcall(m, "V", "addPatientsFromBuffer", paste(match_header, match_data, sep = "\n"))

    ## Run
    tryCatch(.jcall(obj = m, returnSig = "V", method = "run"),
             NumberFormatException = function(e) {
                 e$jobj$printStackTrace()
             })

    ## Load data into R
    matches <- read.table(outfilePath, header=TRUE, sep="\t")
    matched_cohort <- merge(x    = matches,
                            y    = cohort,
                            by.x = "pat_id",
                            by.y = names(cohort)[1], # ID
                            all  = FALSE)

    ## Java null object reference
    m <- .jnull("Object")

    ## Return
    matched_cohort[with(matched_cohort, order(set_num, group_indicator)), ]
}



### Caliper determining function
CalcCaliper3Way <- function(ps1, ps2, group) {

    caliper3Way <- cbind(tapply(ps1, group, var),
                         tapply(ps2, group, var)) %>%
    rowMeans %>%
    sum %>%
    `/`(., 3) %>%
    sqrt %>%
    `*`(., 0.6)

    caliper3Way
}


### Safe matcher ensuring the smallest trio count
## data is the full dataset, cohort is the minimal variables
## Expects a specific data structure
Match3WayAugmentData <- function(data, caliper, minTrioCount = 2, caliperFactor = 1.5,
                                 idVar = "ID", txVar = "Tr", twoPs = c("PS_0","PS_1"),
                                 fun = Match3WayRassen) {

    ## Initialize loop variables
    matchedN     <- 0
    caliperPower <- 0 # Starts with caliperFactor^0

    ## Keep looping while we do not have enough sample size
    ## Sample size is 3 times the trio count
    while (matchedN < minTrioCount * 3) {
        ## Perform three-way matching
        matchingInfo <- fun(cohort = data[c(idVar, txVar, twoPs)],
                            caliper = caliper * caliperFactor^caliperPower,
                            numGroups = 3, groupNames = c(0,1,2))

        ## Update loop variables
        matchedN     <- nrow(matchingInfo)
        caliperPower <- caliperPower + 1
    }

    ## When conditions are met, create an augmented dataset and return it
    outData <- merge(x = data,
                     y = matchingInfo[1:2],
                     by.x = idVar, by.y = "pat_id",
                     all.x = TRUE,  all.y = FALSE)

    ## Return as a list containing additional information
    list(outData      = outData,
         caliperPower = caliperPower - 1,
         caliper      = caliper * caliperFactor^(caliperPower - 1))
}


###
### rJava memo
################################################################################

### Look for available methods for version 2.4.15
## .jmethods(m)
## https://www.rforge.net/doc/packages/rJava/jreflection.html
##
##  [1] "public void org.drugepi.match.Match.run() throws java.lang.Exception"
##  [2] "public void org.drugepi.match.Match.initMatch(java.lang.String,int) throws java.lang.Exception"
##  [3] "public void org.drugepi.match.Match.initMatch(org.drugepi.match.Match$MatchType,int) throws java.lang.Exception"
##  [4] "public void org.drugepi.match.Match.initMatch(java.lang.String) throws java.lang.Exception"
##  [5] "public void org.drugepi.match.Match.initMatch(java.lang.String,java.lang.String) throws java.lang.Exception"
##  [6] "public void org.drugepi.match.Match.addMatchGroup(java.lang.String)"
##  [7] "public void org.drugepi.match.Match.addPatients(org.drugepi.util.RowReader) throws java.lang.Exception"
##  [8] "public java.lang.String org.drugepi.match.Match.getOutfilePath()"
##  [9] "public void org.drugepi.match.Match.setOutfilePath(java.lang.String)"
## [10] "public org.drugepi.match.Match$MatchType org.drugepi.match.Match.getMatchType()"
## [11] "public java.lang.String org.drugepi.match.Match.getMatchOutputData()"
## [12] "public void org.drugepi.PharmacoepiTool.addPatients(java.lang.String) throws java.lang.Exception"
## [13] "public void org.drugepi.PharmacoepiTool.addPatients(java.lang.String,java.lang.String,java.lang.String,java.lang.String,java.lang.String) throws java.lang.Exception"
## [14] "public java.lang.String org.drugepi.PharmacoepiTool.getVersion()"
## [15] "public void org.drugepi.PharmacoepiTool.addPatientsFromBuffer(java.lang.String) throws java.lang.Exception"
## [16] "public java.lang.String org.drugepi.PharmacoepiTool.getDescription()"
## [17] "public final void java.lang.Object.wait(long,int) throws java.lang.InterruptedException"
## [18] "public final native void java.lang.Object.wait(long) throws java.lang.InterruptedException"
## [19] "public final void java.lang.Object.wait() throws java.lang.InterruptedException"
## [20] "public boolean java.lang.Object.equals(java.lang.Object)"
## [21] "public java.lang.String java.lang.Object.toString()"
## [22] "public native int java.lang.Object.hashCode()"
## [23] "public final native java.lang.Class java.lang.Object.getClass()"
## [24] "public final native void java.lang.Object.notify()"
## [25] "public final native void java.lang.Object.notifyAll()"


### Look for fields (data storage)
## .jfields(m)
## https://www.rforge.net/doc/packages/rJava/jreflection.html
##
##  [1] "public final java.lang.String org.drugepi.match.Match.description"
##  [2] "public java.lang.String org.drugepi.match.Match.outfilePath"
##  [3] "public static final double org.drugepi.match.Match.INVALID_CALIPER"
##  [4] "public int org.drugepi.match.Match.matchRatio"
##  [5] "public int org.drugepi.match.Match.fixedRatio"
##  [6] "public double org.drugepi.match.Match.caliper"
##  [7] "public int org.drugepi.match.Match.parallelMatchingMode"
##  [8] "public int org.drugepi.match.Match.startDigit"
##  [9] "public int org.drugepi.match.Match.endDigit"
## [10] "public static java.lang.String org.drugepi.PharmacoepiTool.description"
## [11] "public static java.lang.String org.drugepi.PharmacoepiTool.version"


### Exception handling
## .jclear(), .jcheck()
## https://www.rforge.net/doc/packages/rJava/html/jcheck.html


### Nullify object
## .jnull(m)
## https://www.rforge.net/doc/packages/rJava/html/jnull.html
