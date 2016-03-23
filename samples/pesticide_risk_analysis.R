#:-----------------------------------------------------------------------------:
# Pesticide Risk Analysis Example
#
# NOTE: This is a work in progress. It is in need of refactoring to make it
#       more generic and reusable. Currently, there is a high degree of 
#       reliance on hardcoded values and variables relating to the specific
#       datasets used. The intent is to rewrite this generically to work with 
#       new datasets. Further, during analysis, other tools were used such as
#       IBM Crystal Ball and EPA BMDS. Future vesions of this program will
#       not depend upon the use of external tools, such as these, but would
#       only use R code contained in the script, plus R packages called by it.
# 
# MIT License
#
# Copyright (c) 2016 Jane G. Pouzou and Brian High
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#:-----------------------------------------------------------------------------:

#:-----------------------------------------------------------------------------:
# Clear workspace and load packages
#:-----------------------------------------------------------------------------:

# Clear the workspace, unless you are running in knitr context.
if (!isTRUE(getOption('knitr.in.progress'))) {
    closeAllConnections()
    rm(list=ls())
}

# Load one or more packages into memory, installing as needed.
load.pkgs <- function(pkgs, repos="http://cran.r-project.org") {
    result <- sapply(pkgs, function(pkg) { 
        if (!suppressWarnings(require(pkg, character.only=TRUE))) {
            install.packages(pkg, quiet=TRUE, repos=repos)
            library(pkg, character.only=TRUE)}})
}

# Load the required packages, installing as needed.
load.pkgs(c("dplyr", "mc2d", "fitdistrplus", "STAND", "ggplot2"))

#:-----------------------------------------------------------------------------:
# Import the data
#:-----------------------------------------------------------------------------:

# Import of exposure datasets.
ocap <- read.delim(file.path('data', 'ocapd.txt'))
dfml <- read.delim(file.path('data', 'dfmld.txt'))

# Count numbers of observations per study. See: "Distribution weights", below.
count(dfml, "study")   # 5 studies, with 5 observations per study
count(ocair, "study")  # 4 studies, with 15, 3, 5, and 5 observations

#:-----------------------------------------------------------------------------:
# Normalize exposures by dividing by weight of active ingredient handled
#:-----------------------------------------------------------------------------:

dfml.norm <- data.frame(sapply(names(dfml)[5:12], 
                               function(x) dfml[[x]]/dfml$ailbs))
names(dfml.norm) <- c('normlowa', 'normupa', 'normchest', 'normback', 
                      'normlowl', 'normupl', 'normhead', 'normhand')
dfml <- cbind(dfml, dfml.norm)

ocap.norm <- data.frame(sapply(names(ocap)[c(4:12, 14)], 
                               function(x) ocap[[x]]/ocap$ailbs))
names(ocap.norm) <- c('normlowa', 'normupa', 'normchest', 'normback', 
                      'normlowl', 'normupl', 'normface', 'normhand', 
                      'normCRhead', 'normhead')
ocap <- cbind(ocap, ocap.norm)

#:-----------------------------------------------------------------------------:
# Fit Inhalation Distributions
#:-----------------------------------------------------------------------------:

# (non - nested)
dfair <- subset(dfml, !is.na(dfml$airsamp))
ocair <- subset(ocap, !is.na(ocap$airsamp))
airdistdf <- fitdist(dfair$airsamp, "lnorm")
airdistoc <- fitdist(ocair$airsamp, "lnorm")

# Mixing and loading of Dry Flowable (DF) pesticide task
airdist17 <- fitdist(dfair[dfair$study=="AHE17", "airsamp"], "lnorm")
airdist18 <- fitdist(dfair[dfair$study=="AHE18", "airsamp"], "lnorm")
airdist20 <- fitdist(dfair[dfair$study=="AHE20", "airsamp"], "lnorm")
airdist21 <- fitdist(dfair[dfair$study=="AHE21", "airsamp"], "lnorm")

# Open Cab (OC) application task
airdist07 <- fitdist(ocair[ocair$study=="AHE07", "airsamp"], "lnorm")
airdist62 <- fitdist(ocair[ocair$study=="AHE62", "airsamp"], "lnorm")
airdist63 <- fitdist(ocair[ocair$study=="AHE63", "airsamp"], "lnorm")
airdist64 <- fitdist(ocair[ocair$study=="AHE64", "airsamp"], "lnorm")

# Alternate: Store in dataframes containing lists
# This way, you do not have to hardcode study IDs in variable names.
# And this makes the code more generic, and therefore more reusable.
airfitdf <- data.frame(t(sapply(unique(dfair$study), function(x) 
    fitdist(dfair[dfair$study==x, "airsamp"], "lnorm"))))
airfitdf$study <- unique(dfair$study)
# To access "meanlog" for "AHE17":
# unlist(airfitdf[airfitdf$study=="AHE17", "estimate"])["meanlog"]

airfitoc <- data.frame(t(sapply(unique(ocair$study), function(x) 
    fitdist(ocair[ocair$study==x, "airsamp"], "lnorm"))))
airfitoc$study <- unique(ocair$study)

#:-----------------------------------------------------------------------------:
# Fit Dermal Distributions
#:-----------------------------------------------------------------------------:

# Dry flowables
dflowl <- subset(dfml, !is.na(dfml$lowlegug))
dfupl <- subset(dfml, !is.na(dfml$uplegug))
dflowa <- subset(dfml, !is.na(dfml$lowarmug))
dfupa <- subset(dfml, !is.na(dfml$uparmug))
dfchest <- subset(dfml, !is.na(dfml$chestug))
dfback <- subset(dfml, !is.na(dfml$backug))
dfhand <- subset(dfml, !is.na(dfml$handug))
dfhead <- subset(dfml, !is.na(dfml$faceug))

lowldistdf <- fitdist(dflowl$normlowl, "lnorm")
upldistdf <- fitdist(dfupl$normupl, "lnorm")
upadistdf <- fitdist(dfupa$normupa, "lnorm")
lowadistdf <- fitdist(dflowa$normlowa, "lnorm")
chestdistdf <- fitdist(dfchest$normchest, "lnorm")
backdistdf <- fitdist(dfback$normback, "lnorm")
handdistdf <- fitdist(dfhand$normhand, "lnorm")
headdistdf <- fitdist(dfhead$normhead, "lnorm")

lowldist17 <- fitdist(dflowl[dflowl$study=="AHE17", "normlowl"], "lnorm")
lowldist18 <- fitdist(dflowl[dflowl$study=="AHE18", "normlowl"], "lnorm")
lowldist20 <- fitdist(dflowl[dflowl$study=="AHE20", "normlowl"], "lnorm")
lowldist21 <- fitdist(dflowl[dflowl$study=="AHE21", "normlowl"], "lnorm")

upldist17 <- fitdist(dfupl[dfupl$study=="AHE17", "normupl"], "lnorm")
upldist18 <- fitdist(dfupl[dfupl$study=="AHE18", "normupl"], "lnorm")
upldist20 <- fitdist(dfupl[dfupl$study=="AHE20", "normupl"], "lnorm")
upldist21 <- fitdist(dfupl[dfupl$study=="AHE21", "normupl"], "lnorm")

upadist17 <- fitdist(dfupa[dfupa$study=="AHE17", "normupa"], "lnorm")
upadist18 <- fitdist(dfupa[dfupa$study=="AHE18", "normupa"], "lnorm")
upadist20 <- fitdist(dfupa[dfupa$study=="AHE20", "normupa"], "lnorm")
upadist21 <- fitdist(dfupa[dfupa$study=="AHE21", "normupa"], "lnorm")

lowadist17 <- fitdist(dflowa[dflowa$study=="AHE17", "normlowa"], "lnorm")
lowadist18 <- fitdist(dflowa[dflowa$study=="AHE18", "normlowa"], "lnorm")
lowadist20 <- fitdist(dflowa[dflowa$study=="AHE20", "normlowa"], "lnorm")
lowadist21 <- fitdist(dflowa[dflowa$study=="AHE21", "normlowa"], "lnorm")

chestdist17 <- fitdist(dfchest[dfchest$study=="AHE17", "normchest"], "lnorm")
chestdist18 <- fitdist(dfchest[dfchest$study=="AHE18", "normchest"], "lnorm")
chestdist20 <- fitdist(dfchest[dfchest$study=="AHE20", "normchest"], "lnorm")
chestdist21 <- fitdist(dfchest[dfchest$study=="AHE21", "normchest"], "lnorm")

backdist17 <- fitdist(dfback[dfback$study=="AHE17", "normback"], "lnorm")
backdist18 <- fitdist(dfback[dfback$study=="AHE18", "normback"], "lnorm")
backdist20 <- fitdist(dfback[dfback$study=="AHE20", "normback"], "lnorm")
backdist21 <- fitdist(dfback[dfback$study=="AHE21", "normback"], "lnorm")

handdist17 <- fitdist(dfhand[dfhand$study=="AHE17", "normhand"], "lnorm")
handdist18 <- fitdist(dfhand[dfhand$study=="AHE18", "normhand"], "lnorm")
handdist20 <- fitdist(dfhand[dfhand$study=="AHE20", "normhand"], "lnorm")
handdist21 <- fitdist(dfhand[dfhand$study=="AHE21", "normhand"], "lnorm")

headdist17 <- fitdist(dfhead[dfhead$study=="AHE17", "normhead"], "lnorm")
headdist18 <- fitdist(dfhead[dfhead$study=="AHE18", "normhead"], "lnorm")
headdist20 <- fitdist(dfhead[dfhead$study=="AHE20", "normhead"], "lnorm")
headdist21 <- fitdist(dfhead[dfhead$study=="AHE21", "normhead"], "lnorm")

# Open cab appl
ocair <- subset(ocap, !is.na(ocap$airsamp))
oclowl <- subset(ocap, !is.na(ocap$lowlegug))
ocupl <- subset(ocap, !is.na(ocap$uplegug))
oclowa <- subset(ocap, !is.na(ocap$lowarmug))
ocupa <- subset(ocap, !is.na(ocap$uparmug))
occhest <- subset(ocap, !is.na(ocap$chestug))
ocback <- subset(ocap, !is.na(ocap$backug))
ochand <- subset(ocap, !is.na(ocap$handug))
ochead <- subset(ocap, !is.na(ocap$headug))
ocface <- subset(ocap, !is.na(ocap$faceug))
ocCRhead <- subset(ocap, !is.na(ocap$CRheadug))

lowl07 <- subset(oclowl, oclowl$study=="AHE07")
lowl62 <- subset(oclowl, oclowl$study=="AHE62")
lowl63 <- subset(oclowl, oclowl$study=="AHE63")
lowl64 <- subset(oclowl, oclowl$study=="AHE64")

lowldistoc <- fitdist(oclowl$normlowl, "lnorm")
lowldist07 <- fitdist(lowl07$normlowl, "lnorm")
lowldist62 <- fitdist(lowl62$normlowl, "lnorm")
lowldist63 <- fitdist(lowl63$normlowl, "lnorm")
lowldist64 <- fitdist(lowl64$normlowl, "lnorm")

upl07 <- subset(ocupl, ocupl$study=="AHE07")
upl62 <- subset(ocupl, ocupl$study=="AHE62")
upl63 <- subset(ocupl, ocupl$study=="AHE63")
upl64 <- subset(ocupl, ocupl$study=="AHE64")

upldistoc <- fitdist(ocupl$normupl, "lnorm")
upldist07 <- fitdist(upl07$normupl, "lnorm")
upldist62 <- fitdist(upl62$normupl, "lnorm")
upldist63 <- fitdist(upl63$normupl, "lnorm")
upldist64 <- fitdist(upl64$normupl, "lnorm")

upa07 <- subset(ocupa, ocupa$study=="AHE07")
upa62 <- subset(ocupa, ocupa$study=="AHE62")
upa63 <- subset(ocupa, ocupa$study=="AHE63")
upa64 <- subset(ocupa, ocupa$study=="AHE64")

upadistoc <- fitdist(ocupa$normupa, "lnorm")
upadist07 <- fitdist(upa07$normupa, "lnorm")
upadist62 <- fitdist(upa62$normupa, "lnorm")
upadist63 <- fitdist(upa63$normupa, "lnorm")
upadist64 <- fitdist(upa64$normupa, "lnorm")

lowa07 <- subset(oclowa, oclowa$study=="AHE07")
lowa62 <- subset(oclowa, oclowa$study=="AHE62")
lowa63 <- subset(oclowa, oclowa$study=="AHE63")
lowa64 <- subset(oclowa, oclowa$study=="AHE64")

lowadistoc <- fitdist(oclowa$normlowa, "lnorm")
lowadist07 <- fitdist(lowa07$normlowa, "lnorm")
lowadist62 <- fitdist(lowa62$normlowa, "lnorm")
lowadist63 <- fitdist(lowa63$normlowa, "lnorm")
lowadist64 <- fitdist(lowa64$normlowa, "lnorm")

chest07 <- subset(occhest, occhest$study=="AHE07")
chest62 <- subset(occhest, occhest$study=="AHE62")
chest63 <- subset(occhest, occhest$study=="AHE63")
chest64 <- subset(occhest, occhest$study=="AHE64")

chestdistoc <- fitdist(occhest$normchest, "lnorm")
chestdist07 <- fitdist(chest07$normchest, "lnorm")
chestdist62 <- fitdist(chest62$normchest, "lnorm")
chestdist63 <- fitdist(chest63$normchest, "lnorm")
chestdist64 <- fitdist(chest64$normchest, "lnorm")

back07 <- subset(ocback, ocback$study=="AHE07")
back62 <- subset(ocback, ocback$study=="AHE62")
back63 <- subset(ocback, ocback$study=="AHE63")
back64 <- subset(ocback, ocback$study=="AHE64")

backdistoc <- fitdist(ocback$normback, "lnorm")
backdist07 <- fitdist(back07$normback, "lnorm")
backdist62 <- fitdist(back62$normback, "lnorm")
backdist63 <- fitdist(back63$normback, "lnorm")
backdist64 <- fitdist(back64$normback, "lnorm")

hand07 <- subset(ochand, ochand$study=="AHE07")
hand62 <- subset(ochand, ochand$study=="AHE62")
hand63 <- subset(ochand, ochand$study=="AHE63")
hand64 <- subset(ochand, ochand$study=="AHE64")

handdistoc <- fitdist(ochand$normhand, "lnorm")
handdist07 <- fitdist(hand07$normhand, "lnorm")
handdist62 <- fitdist(hand62$normhand, "lnorm")
handdist63 <- fitdist(hand63$normhand, "lnorm")
handdist64 <- fitdist(hand64$normhand, "lnorm")

head07 <- subset(ochead, ochead$study=="AHE07")
head62 <- subset(ochead, ochead$study=="AHE62")
head63 <- subset(ochead, ochead$study=="AHE63")
head64 <- subset(ochead, ochead$study=="AHE64")

headdistoc <- fitdist(ochead$normhead, "lnorm")
headdist07 <- fitdist(head07$normhead, "lnorm")
headdist62 <- fitdist(head62$normhead, "lnorm")
headdist63 <- fitdist(head63$normhead, "lnorm")
headdist64 <- fitdist(head64$normhead, "lnorm")

CRhead07 <- subset(ocCRhead, ocCRhead$study=="AHE07")
CRhead62 <- subset(ocCRhead, ocCRhead$study=="AHE62")
CRhead63 <- subset(ocCRhead, ocCRhead$study=="AHE63")
CRhead64 <- subset(ocCRhead, ocCRhead$study=="AHE64")

CRheaddistoc <- fitdist(ocCRhead$normCRhead, "lnorm")
CRheaddist07 <- fitdist(CRhead07$normCRhead, "lnorm")
CRheaddist62 <- fitdist(CRhead62$normCRhead, "lnorm")
CRheaddist63 <- fitdist(CRhead63$normCRhead, "lnorm")
CRheaddist64 <- fitdist(CRhead64$normCRhead, "lnorm")

face07 <- subset(ocface, ocface$study=="AHE07")
face62 <- subset(ocface, ocface$study=="AHE62")
face63 <- subset(ocface, ocface$study=="AHE63")
face64 <- subset(ocface, ocface$study=="AHE64")

facedistoc <- fitdist(ocface$normface, "lnorm")
facedist07 <- fitdist(face07$normface, "lnorm")
facedist62 <- fitdist(face62$normface, "lnorm")
facedist63 <- fitdist(face63$normface, "lnorm")
facedist64 <- fitdist(face64$normface, "lnorm")

ndvar=(1001)
ndunc=(101)

#:-----------------------------------------------------------------------------:
# Define distributions
#:-----------------------------------------------------------------------------:

set.seed(1)

# Acres of application (Source, AHETF source doc)
acredis <- mcstoc(rtriang, type="V", min=1, max=60, mode=40)

# Application Error distribution (source:  Rider and Dickey, 1982)
apprerr <- mcstoc(rnorm, type="U", mean=-0.028, sd=0.262)

# Body weight in kilograms (source: EFH, taken from NHANES IV, males 18-65, 
# and from participants in AHETF studies).
bwdis <- mcstoc(rnorm, type="V", mean=85.47, sd=19.03, rtrunc=TRUE, linf=0, 
                lsup=300)
bwdis2 <- mcstoc(rnorm, type="V", mean=87.25, sd=16.84, rtrunc=TRUE, linf=0, 
                 lsup=300)
bdwgt <- mcstoc(rempiricalD, type="V", values=1018:1019, prob=c(0.5, 0.5))

#:-----------------------------------------------------------------------------:
# Breathing Rate (Source:  Exposure Factors Handbook)
#:-----------------------------------------------------------------------------:

# Applicator Breathing Rate in m^3/hr (Sedentary)
applbr1 <- mcstoc(rnorm, type="V", mean=.071, sd=0.4, rtrunc=TRUE, linf=0)
applbr2 <- mcstoc(rnorm, type="V", mean=.078, sd=0.36, rtrunc=TRUE, linf=0)
applbr3 <- mcstoc(runif, type="V", min=0.5, max=0.5)
brwgt <- mcstoc(rempiricalD, type="V", values=1020:1022, prob=c(1/3, 1/3, 1/3))

# Mixer/Loader Breathing Rate m^3/hr (Moderate Activity)
mlbr1 <- mcstoc(rnorm, type="V", mean=.84, sd=0.47, rtrunc=TRUE, linf=0)
mlbr2 <- mcstoc(rnorm, type="V", mean=.84, sd=0.54, rtrunc=TRUE, linf=0)
mlbr3 <- mcstoc(runif, type="V", min=1, max=1)

# Skin surface area by body part (Exposure Factors handbook, for males age 21+)
lowlegsa <- mcstoc(rnorm, type="V", mean=2680, sd=340.48, rtrunc=TRUE, linf=0)
uplegsa <- mcstoc(rnorm, type="V", mean=4120, sd=674.87, rtrunc=TRUE, linf=0)
chestsa <- mcstoc(rnorm, type="V", mean=3875, sd=829.90, rtrunc=TRUE, linf=0)
backsa <- mcstoc(rnorm, type="V", mean=3875, sd=829.90, rtrunc=TRUE, linf=0)
uparmsa <- mcstoc(rnorm, type="V", mean=1720, sd=291.84, rtrunc=TRUE, linf=0)
lowarmsa <- mcstoc(rnorm, type="V", mean=1480, sd=297.92, rtrunc=TRUE, linf=0)
headnecksa <- mcstoc(rnorm, type="V", mean=1620, sd=109.44, rtrunc=TRUE, linf=0)
facefnecksa <- mcstoc(rnorm, type="V", mean=583, sd=36.48, rtrunc=TRUE, linf=0)

#:-----------------------------------------------------------------------------:
# Distribution weights
#:-----------------------------------------------------------------------------:

# DF weights:  n=5 for each study
dfweight <- mcstoc(rempiricalD, type="V", values=1:4, 
                   prob=c(0.25, 0.25, 0.25, 0.25))

# OC weights
# 07, n=15; 62, n=3; 63, n=5;  64, n=5
#
# totN <- 15+3+5+5
# 15/totN
## 0.5357143
# 3/totN
## 0.1071429
# 5/totN
## 0.1785714
#
ocwgt <- mcstoc(rempiricalD, type="V", values=87:90, 
                prob=c(0.5357, 0.107, 0.17, 0.17))

#:-----------------------------------------------------------------------------:
# Inhalation study distributions as mcnodes
#:-----------------------------------------------------------------------------:

# oc nested
mcairdist07 <- mcstoc(rlnorm, type="V", meanlog=airdist07$est["meanlog"], 
                      sdlog=airdist07$est["sdlog"], rtrunc=TRUE, linf=0)
mcairdist62 <- mcstoc(rlnorm, type="V", meanlog=airdist62$est["meanlog"], 
                      sdlog=airdist62$est["sdlog"], rtrunc=TRUE, linf=0) 
mcairdist63 <- mcstoc(rlnorm, type="V", meanlog=airdist63$est["meanlog"], 
                      sdlog=airdist63$est["sdlog"], rtrunc=TRUE, linf=0)  
mcairdist64 <- mcstoc(rlnorm, type="V", meanlog=airdist64$est["meanlog"], 
                      sdlog=airdist64$est["sdlog"], rtrunc=TRUE, linf=0)

# df nested
mcairdist17 <- mcstoc(rlnorm, type="V", meanlog=airdist17$est["meanlog"], 
                      sdlog=airdist17$est["sdlog"], rtrunc=TRUE, linf=0)
mcairdist18 <- mcstoc(rlnorm, type="V", meanlog=airdist18$est["meanlog"], 
                      sdlog=airdist18$est["sdlog"], rtrunc=TRUE, linf=0)  
mcairdist20 <- mcstoc(rlnorm, type="V", meanlog=airdist20$est["meanlog"], 
                      sdlog=airdist20$est["sdlog"], rtrunc=TRUE, linf=0) 
mcairdist21 <- mcstoc(rlnorm, type="V", meanlog=airdist21$est["meanlog"], 
                      sdlog=airdist21$est["sdlog"], rtrunc=TRUE, linf=0)

#:-----------------------------------------------------------------------------:
# Dermal study distributions as mcnodes 
#:-----------------------------------------------------------------------------:

# Open Cab Appl
mclowldist07 <- mcstoc(rlnorm, type="V", meanlog=lowldist07$est["meanlog"], 
                       sdlog=lowldist07$est["sdlog"], rtrunc=TRUE, linf=0)
mclowldist62 <- mcstoc(rlnorm, type="V", meanlog=lowldist62$est["meanlog"], 
                       sdlog=lowldist62$est["sdlog"], rtrunc=TRUE, linf=0)
mclowldist63 <- mcstoc(rlnorm, type="V", meanlog=lowldist63$est["meanlog"], 
                       sdlog=lowldist63$est["sdlog"], rtrunc=TRUE, linf=0)
mclowldist64 <- mcstoc(rlnorm, type="V", meanlog=lowldist64$est["meanlog"], 
                       sdlog=lowldist64$est["sdlog"], rtrunc=TRUE, linf=0)
mcupldist07 <- mcstoc(rlnorm, type="V", meanlog=upldist07$est["meanlog"], 
                      sdlog=upldist07$est["sdlog"], rtrunc=TRUE, linf=0)
mcupldist62 <- mcstoc(rlnorm, type="V", meanlog=upldist62$est["meanlog"], 
                      sdlog=upldist62$est["sdlog"], rtrunc=TRUE, linf=0)
mcupldist63 <- mcstoc(rlnorm, type="V", meanlog=upldist63$est["meanlog"], 
                      sdlog=upldist63$est["sdlog"], rtrunc=TRUE, linf=0)
mcupldist64 <- mcstoc(rlnorm, type="V", meanlog=upldist64$est["meanlog"], 
                      sdlog=upldist64$est["sdlog"], rtrunc=TRUE, linf=0)
mclowadist07 <- mcstoc(rlnorm, type="V", meanlog=lowadist07$est["meanlog"], 
                       sdlog=lowadist07$est["sdlog"], rtrunc=TRUE, linf=0)
mclowadist62 <- mcstoc(rlnorm, type="V", meanlog=lowadist62$est["meanlog"], 
                       sdlog=lowadist62$est["sdlog"], rtrunc=TRUE, linf=0)
mclowadist63 <- mcstoc(rlnorm, type="V", meanlog=lowadist63$est["meanlog"], 
                       sdlog=lowadist63$est["sdlog"], rtrunc=TRUE, linf=0)
mclowadist64 <- mcstoc(rlnorm, type="V", meanlog=lowadist64$est["meanlog"], 
                       sdlog=lowadist64$est["sdlog"], rtrunc=TRUE, linf=0)
mcupadist07 <- mcstoc(rlnorm, type="V", meanlog=upadist07$est["meanlog"], 
                      sdlog=upadist07$est["sdlog"], rtrunc=TRUE, linf=0)
mcupadist62 <- mcstoc(rlnorm, type="V", meanlog=upadist62$est["meanlog"], 
                      sdlog=upadist62$est["sdlog"], rtrunc=TRUE, linf=0)
mcupadist63 <- mcstoc(rlnorm, type="V", meanlog=upadist63$est["meanlog"], 
                      sdlog=upadist63$est["sdlog"], rtrunc=TRUE, linf=0)
mcupadist64 <- mcstoc(rlnorm, type="V", meanlog=upadist64$est["meanlog"], 
                      sdlog=upadist64$est["sdlog"], rtrunc=TRUE, linf=0)
mcchestdist07 <- mcstoc(rlnorm, type="V", meanlog=chestdist07$est["meanlog"], 
                        sdlog=chestdist07$est["sdlog"], rtrunc=TRUE, linf=0)
mcchestdist62 <- mcstoc(rlnorm, type="V", meanlog=chestdist62$est["meanlog"], 
                        sdlog=chestdist62$est["sdlog"], rtrunc=TRUE, linf=0)
mcchestdist63 <- mcstoc(rlnorm, type="V", meanlog=chestdist63$est["meanlog"], 
                        sdlog=chestdist63$est["sdlog"], rtrunc=TRUE, linf=0)
mcchestdist64 <- mcstoc(rlnorm, type="V", meanlog=chestdist64$est["meanlog"], 
                        sdlog=chestdist64$est["sdlog"], rtrunc=TRUE, linf=0)
mcbackdist07 <- mcstoc(rlnorm, type="V", meanlog=backdist07$est["meanlog"], 
                       sdlog=backdist07$est["sdlog"], rtrunc=TRUE, linf=0)
mcbackdist62 <- mcstoc(rlnorm, type="V", meanlog=backdist62$est["meanlog"], 
                       sdlog=backdist62$est["sdlog"], rtrunc=TRUE, linf=0)
mcbackdist63 <- mcstoc(rlnorm, type="V", meanlog=backdist63$est["meanlog"], 
                       sdlog=backdist63$est["sdlog"], rtrunc=TRUE, linf=0)
mcbackdist64 <- mcstoc(rlnorm, type="V", meanlog=backdist64$est["meanlog"], 
                       sdlog=backdist64$est["sdlog"], rtrunc=TRUE, linf=0)
mcheaddist07 <- mcstoc(rlnorm, type="V", meanlog=headdist07$est["meanlog"], 
                       sdlog=headdist07$est["sdlog"], rtrunc=TRUE, linf=0)
mcheaddist62 <- mcstoc(rlnorm, type="V", meanlog=headdist62$est["meanlog"], 
                       sdlog=headdist62$est["sdlog"], rtrunc=TRUE, linf=0)
mcheaddist63 <- mcstoc(rlnorm, type="V", meanlog=headdist63$est["meanlog"], 
                       sdlog=headdist63$est["sdlog"], rtrunc=TRUE, linf=0)
mcheaddist64 <- mcstoc(rlnorm, type="V", meanlog=headdist64$est["meanlog"], 
                       sdlog=headdist64$est["sdlog"], rtrunc=TRUE, linf=0)
mchanddist07 <- mcstoc(rlnorm, type="V", meanlog=handdist07$est["meanlog"], 
                       sdlog=handdist07$est["sdlog"], rtrunc=TRUE, linf=0)
mchanddist62 <- mcstoc(rlnorm, type="V", meanlog=handdist62$est["meanlog"], 
                       sdlog=handdist62$est["sdlog"], rtrunc=TRUE, linf=0)
mchanddist63 <- mcstoc(rlnorm, type="V", meanlog=handdist63$est["meanlog"], 
                       sdlog=handdist63$est["sdlog"], rtrunc=TRUE, linf=0)
mchanddist64 <- mcstoc(rlnorm, type="V", meanlog=handdist64$est["meanlog"], 
                       sdlog=handdist64$est["sdlog"], rtrunc=TRUE, linf=0)
mcfacedist07 <- mcstoc(rlnorm, type="V", meanlog=facedist07$est["meanlog"], 
                       sdlog=facedist07$est["sdlog"], rtrunc=TRUE, linf=0)
mcfacedist62 <- mcstoc(rlnorm, type="V", meanlog=facedist62$est["meanlog"], 
                       sdlog=facedist62$est["sdlog"], rtrunc=TRUE, linf=0)
mcfacedist63 <- mcstoc(rlnorm, type="V", meanlog=facedist63$est["meanlog"], 
                       sdlog=facedist63$est["sdlog"], rtrunc=TRUE, linf=0)
mcfacedist64 <- mcstoc(rlnorm, type="V", meanlog=facedist64$est["meanlog"], 
                       sdlog=facedist64$est["sdlog"], rtrunc=TRUE, linf=0)
mcCRheaddist07 <- mcstoc(rlnorm, type="V", meanlog=CRheaddist07$est["meanlog"], 
                         sdlog=CRheaddist07$est["sdlog"], rtrunc=TRUE, linf=0)
mcCRheaddist62 <- mcstoc(rlnorm, type="V", meanlog=CRheaddist62$est["meanlog"], 
                         sdlog=CRheaddist62$est["sdlog"], rtrunc=TRUE, linf=0)
mcCRheaddist63 <- mcstoc(rlnorm, type="V", meanlog=CRheaddist63$est["meanlog"], 
                         sdlog=CRheaddist63$est["sdlog"], rtrunc=TRUE, linf=0)
mcCRheaddist64 <- mcstoc(rlnorm, type="V", meanlog=CRheaddist64$est["meanlog"], 
                         sdlog=CRheaddist64$est["sdlog"], rtrunc=TRUE, linf=0)

# Dry flowables
mclowldist17 <- mcstoc(rlnorm, type="V", meanlog=lowldist17$est["meanlog"], 
                       sdlog=lowldist17$est["sdlog"], rtrunc=TRUE, linf=0)
mclowldist18 <- mcstoc(rlnorm, type="V", meanlog=lowldist18$est["meanlog"], 
                       sdlog=lowldist18$est["sdlog"], rtrunc=TRUE, linf=0)
mclowldist20 <- mcstoc(rlnorm, type="V", meanlog=lowldist20$est["meanlog"], 
                       sdlog=lowldist20$est["sdlog"], rtrunc=TRUE, linf=0)
mclowldist21 <- mcstoc(rlnorm, type="V", meanlog=lowldist21$est["meanlog"], 
                       sdlog=lowldist21$est["sdlog"], rtrunc=TRUE, linf=0)
mcupldist17 <- mcstoc(rlnorm, type="V", meanlog=upldist17$est["meanlog"], 
                      sdlog=upldist17$est["sdlog"], rtrunc=TRUE, linf=0)
mcupldist18 <- mcstoc(rlnorm, type="V", meanlog=upldist18$est["meanlog"], 
                      sdlog=upldist18$est["sdlog"], rtrunc=TRUE, linf=0)
mcupldist20 <- mcstoc(rlnorm, type="V", meanlog=upldist20$est["meanlog"], 
                      sdlog=upldist20$est["sdlog"], rtrunc=TRUE, linf=0)
mcupldist21 <- mcstoc(rlnorm, type="V", meanlog=upldist21$est["meanlog"], 
                      sdlog=upldist21$est["sdlog"], rtrunc=TRUE, linf=0)
mclowadist17 <- mcstoc(rlnorm, type="V", meanlog=lowadist17$est["meanlog"], 
                       sdlog=lowadist17$est["sdlog"], rtrunc=TRUE, linf=0)
mclowadist18 <- mcstoc(rlnorm, type="V", meanlog=lowadist18$est["meanlog"], 
                       sdlog=lowadist18$est["sdlog"], rtrunc=TRUE, linf=0)
mclowadist20 <- mcstoc(rlnorm, type="V", meanlog=lowadist20$est["meanlog"], 
                       sdlog=lowadist20$est["sdlog"], rtrunc=TRUE, linf=0)
mclowadist21 <- mcstoc(rlnorm, type="V", meanlog=lowadist21$est["meanlog"], 
                       sdlog=lowadist21$est["sdlog"], rtrunc=TRUE, linf=0)
mcupadist17 <- mcstoc(rlnorm, type="V", meanlog=upadist17$est["meanlog"], 
                      sdlog=upadist17$est["sdlog"], rtrunc=TRUE, linf=0)
mcupadist18 <- mcstoc(rlnorm, type="V", meanlog=upadist18$est["meanlog"], 
                      sdlog=upadist18$est["sdlog"], rtrunc=TRUE, linf=0)
mcupadist20 <- mcstoc(rlnorm, type="V", meanlog=upadist20$est["meanlog"], 
                      sdlog=upadist20$est["sdlog"], rtrunc=TRUE, linf=0)
mcupadist21 <- mcstoc(rlnorm, type="V", meanlog=upadist21$est["meanlog"], 
                      sdlog=upadist21$est["sdlog"], rtrunc=TRUE, linf=0)
mcchestdist17 <- mcstoc(rlnorm, type="V", meanlog=chestdist17$est["meanlog"], 
                        sdlog=chestdist17$est["sdlog"], rtrunc=TRUE, linf=0)
mcchestdist18 <- mcstoc(rlnorm, type="V", meanlog=chestdist18$est["meanlog"], 
                        sdlog=chestdist18$est["sdlog"], rtrunc=TRUE, linf=0)
mcchestdist20 <- mcstoc(rlnorm, type="V", meanlog=chestdist20$est["meanlog"], 
                        sdlog=chestdist20$est["sdlog"], rtrunc=TRUE, linf=0)
mcchestdist21 <- mcstoc(rlnorm, type="V", meanlog=chestdist21$est["meanlog"], 
                        sdlog=chestdist21$est["sdlog"], rtrunc=TRUE, linf=0)
mcbackdist17 <- mcstoc(rlnorm, type="V", meanlog=backdist17$est["meanlog"], 
                       sdlog=backdist17$est["sdlog"], rtrunc=TRUE, linf=0)
mcbackdist18 <- mcstoc(rlnorm, type="V", meanlog=backdist18$est["meanlog"], 
                       sdlog=backdist18$est["sdlog"], rtrunc=TRUE, linf=0)
mcbackdist20 <- mcstoc(rlnorm, type="V", meanlog=backdist20$est["meanlog"], 
                       sdlog=backdist20$est["sdlog"], rtrunc=TRUE, linf=0)
mcbackdist21 <- mcstoc(rlnorm, type="V", meanlog=backdist21$est["meanlog"], 
                       sdlog=backdist21$est["sdlog"], rtrunc=TRUE, linf=0)
mcheaddist17 <- mcstoc(rlnorm, type="V", meanlog=headdist17$est["meanlog"], 
                       sdlog=headdist17$est["sdlog"], rtrunc=TRUE, linf=0)
mcheaddist18 <- mcstoc(rlnorm, type="V", meanlog=headdist18$est["meanlog"], 
                       sdlog=headdist18$est["sdlog"], rtrunc=TRUE, linf=0)
mcheaddist20 <- mcstoc(rlnorm, type="V", meanlog=headdist20$est["meanlog"], 
                       sdlog=headdist20$est["sdlog"], rtrunc=TRUE, linf=0)
mcheaddist21 <- mcstoc(rlnorm, type="V", meanlog=headdist21$est["meanlog"], 
                       sdlog=headdist21$est["sdlog"], rtrunc=TRUE, linf=0)
mchanddist17 <- mcstoc(rlnorm, type="V", meanlog=handdist17$est["meanlog"], 
                       sdlog=handdist17$est["sdlog"], rtrunc=TRUE, linf=0)
mchanddist18 <- mcstoc(rlnorm, type="V", meanlog=handdist18$est["meanlog"], 
                       sdlog=handdist18$est["sdlog"], rtrunc=TRUE, linf=0)
mchanddist20 <- mcstoc(rlnorm, type="V", meanlog=handdist20$est["meanlog"], 
                       sdlog=handdist20$est["sdlog"], rtrunc=TRUE, linf=0)
mchanddist21 <- mcstoc(rlnorm, type="V", meanlog=handdist21$est["meanlog"], 
                       sdlog=handdist21$est["sdlog"], rtrunc=TRUE, linf=0)

#:-----------------------------------------------------------------------------:
# Dermal total nested exposure distributions as mcnodes 
#:-----------------------------------------------------------------------------:

# Open Cab Appl
nmclowldistoc <- mcprobtree(ocwgt, 
                            list("87"=mclowldist07, "88"=mclowldist62, 
                                 "89"=mclowldist63, "90"=mclowldist64)) 
nmcupldistoc <- mcprobtree(ocwgt, 
                           list("87"=mcupldist07, "88"=mcupldist62, 
                                "89"=mcupldist63, "90"=mcupldist64))
nmclowadistoc <- mcprobtree(ocwgt, 
                            list("87"=mclowadist07, "88"=mclowadist62, 
                                 "89"=mclowadist63, "90"=mclowadist64))
nmcupadistoc <- mcprobtree(ocwgt, 
                           list("87"=mcupadist07, "88"=mcupadist62, 
                                "89"=mcupadist63, "90"=mcupadist64))
nmcchestdistoc <- mcprobtree(ocwgt, 
                             list("87"=mcchestdist07, "88"=mcchestdist62, 
                                  "89"=mcchestdist63, "90"=mcchestdist64))
nmcbackdistoc <- mcprobtree(ocwgt, 
                            list("87"=mcbackdist07, "88"=mcbackdist62, 
                                 "89"=mcbackdist63, "90"=mcbackdist64))
nmcheaddistoc <- mcprobtree(ocwgt, 
                            list("87"=mcheaddist07, "88"=mcheaddist62, 
                                 "89"=mcheaddist63, "90"=mcheaddist64))
nmchanddistoc <- mcprobtree(ocwgt, 
                            list("87"=mchanddist07, "88"=mchanddist62, 
                                 "89"=mchanddist63, "90"=mchanddist64))
nmcfacedistoc <- mcprobtree(ocwgt, 
                            list("87"=mcfacedist07, "88"=mcfacedist62, 
                                 "89"=mcfacedist63, "90"=mcfacedist64))
nmcCRheaddistoc <- mcprobtree(ocwgt, 
                              list("87"=mcCRheaddist07, "88"=mcCRheaddist62, 
                                   "89"=mcCRheaddist63, "90"=mcCRheaddist64))

# Dry flowables
nmclowldistdf <- mcprobtree(dfweight, 
                            list("1"=mclowldist17, "2"=mclowldist18, 
                                 "3"=mclowldist20, "4"=mclowldist21))
nmcupldistdf <- mcprobtree(dfweight, 
                           list("1"=mcupldist17, "2"=mcupldist18, 
                                "3"=mcupldist20, "4"=mcupldist21))
nmclowadistdf <- mcprobtree(dfweight, 
                            list("1"=mclowadist17, "2"=mclowadist18, 
                                 "3"=mclowadist20, "4"=mclowadist21))
nmcupadistdf <- mcprobtree(dfweight, 
                           list("1"=mcupadist17, "2"=mcupadist18, 
                                "3"=mcupadist20, "4"=mcupadist21))
nmcchestdistdf <- mcprobtree(dfweight, 
                             list("1"=mcchestdist17, "2"=mcchestdist18, 
                                  "3"=mcchestdist20, "4"=mcchestdist21))
nmcbackdistdf <- mcprobtree(dfweight, 
                            list("1"=mcbackdist17, "2"=mcbackdist18, 
                                 "3"=mcbackdist20, "4"=mcbackdist21))
nmcheaddistdf <- mcprobtree(dfweight, 
                            list("1"=mcheaddist17, "2"=mcheaddist18, 
                                 "3"=mcheaddist20, "4"=mcheaddist21))
nmchanddistdf <- mcprobtree(dfweight, 
                            list("1"=mchanddist17, "2"=mchanddist18, 
                                 "3"=mchanddist20, "4"=mchanddist21))

# Inhalation nested exposure distributions as mcnodes 
bodyweight <- mcprobtree(bdwgt, list("1018"=bwdis, "1019"=bwdis2), type="V")
applbrrate <- mcprobtree(brwgt, type="V", list("1020"=applbr1, "1021"=applbr2, 
                                               "1022"=applbr3))
mlbrrate <- mcprobtree(brwgt, list("1020"=mlbr1, "1021"=mlbr2, "1022"=mlbr3), 
                       type="V")

# oc
nmcairdistoc <- mcprobtree(ocwgt, list("87"=mcairdist07, "88"=mcairdist62, 
                                       "89"=mcairdist63, "90"=mcairdist64))

# df
nmcairdistdf <- mcprobtree(dfweight, list("1"=mcairdist17, "2"=mcairdist18, 
                                          "3"=mcairdist20, "4"=mcairdist21))

#:-----------------------------------------------------------------------------:
# Correlations between body areas' exposure 
#:-----------------------------------------------------------------------------:

# The "target" values came from a separate step performed with IBM Crystal Ball.

# df
cornode(nmclowldistdf, nmcupldistdf, target=0.85)
cornode(nmclowldistdf, nmclowadistdf, target=0.72)
cornode(nmclowldistdf, nmcupadistdf, target=0.63)
cornode(nmclowldistdf, nmcchestdistdf, target=0.76)
cornode(nmclowldistdf, nmcbackdistdf, target=0.63)
cornode(nmclowldistdf, nmcheaddistdf, target=0.56)
cornode(nmclowldistdf, nmchanddistdf, target=0.52)

cornode(nmcupldistdf, nmclowadistdf, target=0.85)
cornode(nmcupldistdf, nmcupadistdf, target=0.74)
cornode(nmcupldistdf, nmcchestdistdf, target=0.90)
cornode(nmcupldistdf, nmcbackdistdf, target=0.73)
cornode(nmcupldistdf, nmcheaddistdf, target=0.65)
cornode(nmcupldistdf, nmchanddistdf, target=0.73)

cornode(nmclowadistdf, nmcupadistdf, target=0.82)
cornode(nmclowadistdf, nmcchestdistdf, target=0.95)
cornode(nmclowadistdf, nmcbackdistdf, target=0.86)
cornode(nmclowadistdf, nmcheaddistdf, target=0.72)
cornode(nmclowadistdf, nmchanddistdf, target=0.58)

cornode(nmcupadistdf, nmcchestdistdf, target=0.82)
cornode(nmcupadistdf, nmcbackdistdf, target=0.91)
cornode(nmcupadistdf, nmcheaddistdf, target=0.75)
cornode(nmcupadistdf, nmchanddistdf, target=0.54)

cornode(nmcchestdistdf, nmcbackdistdf, target=0.82)
cornode(nmcchestdistdf, nmcheaddistdf, target=0.72)
cornode(nmcchestdistdf, nmchanddistdf, target=0.63)

cornode(nmcbackdistdf, nmcheaddistdf, target=0.83)
cornode(nmcbackdistdf, nmchanddistdf, target=0.51)

cornode(nmcheaddistdf, nmchanddistdf, target=0.46)

# Open Cab Appl
cornode(nmclowldistoc, nmclowadistoc, target=0.19)
cornode(nmclowldistoc, nmcupadistoc, target=0.18)
cornode(nmclowldistoc, nmcchestdistoc, target=0.11)
cornode(nmclowldistoc, nmcbackdistoc, target=0.16)
cornode(nmclowldistoc, nmcheaddistoc, target=0.13)
cornode(nmclowldistoc, nmchanddistoc, target=0.26)
cornode(nmclowldistoc, nmcfacedistoc, target=0.12)
cornode(nmclowldistoc, nmcupldistoc, target=0.89)
cornode(nmclowldistoc, nmcCRheaddistoc, target=0.13)

cornode(nmcupldistoc, nmclowadistoc, target=0.20)
cornode(nmcupldistoc, nmcupadistoc, target=0.20)
cornode(nmcupldistoc, nmcchestdistoc, target=0.84)
cornode(nmcupldistoc, nmcbackdistoc, target=0.17)
cornode(nmcupldistoc, nmcheaddistoc, target=0.09)
cornode(nmcupldistoc, nmchanddistoc, target=0.29)
cornode(nmcupldistoc, nmcfacedistoc, target=0.09)
cornode(nmcupldistoc, nmcCRheaddistoc, target=0.09)

cornode(nmclowadistoc, nmcupadistoc, target=0.94)
cornode(nmclowadistoc, nmcchestdistoc, target=0.83)
cornode(nmclowadistoc, nmcbackdistoc, target=0.83)
cornode(nmclowadistoc, nmcheaddistoc, target=0.58)
cornode(nmclowadistoc, nmchanddistoc, target=0.69)
cornode(nmclowadistoc, nmcfacedistoc, target=0.67)
cornode(nmclowadistoc, nmcCRheaddistoc, target=0.58)

cornode(nmcupadistoc, nmcchestdistoc, target=0.83)
cornode(nmcupadistoc, nmcbackdistoc, target=0.8)
cornode(nmcupadistoc, nmcheaddistoc, target=0.59)
cornode(nmcupadistoc, nmchanddistoc, target=0.66)
cornode(nmcupadistoc, nmcfacedistoc, target=0.65)
cornode(nmcupadistoc, nmcCRheaddistoc, target=0.59)

cornode(nmcchestdistoc, nmcbackdistoc, target=0.69)
cornode(nmcchestdistoc, nmcheaddistoc, target=0.7)
cornode(nmcchestdistoc, nmchanddistoc, target=0.65)
cornode(nmcchestdistoc, nmcfacedistoc, target=0.68)
cornode(nmcchestdistoc, nmcCRheaddistoc, target=0.7)

cornode(nmcbackdistoc, nmcheaddistoc, target=0.49)
cornode(nmcbackdistoc, nmchanddistoc, target=0.58)
cornode(nmcbackdistoc, nmcfacedistoc, target=0.8)
cornode(nmcbackdistoc, nmcCRheaddistoc, target=0.49)

cornode(nmcheaddistoc, nmchanddistoc, target=0.46)
cornode(nmcheaddistoc, nmcfacedistoc, target=0.48)
cornode(nmcheaddistoc, nmcCRheaddistoc, target=0.99)

cornode(nmchanddistoc, nmcfacedistoc, target=0.45)
cornode(nmchanddistoc, nmcCRheaddistoc, target=0.46)

cornode(nmcfacedistoc, nmcCRheaddistoc, target=0.48)

#:-----------------------------------------------------------------------------:
# Emamectin benzoate 
#:-----------------------------------------------------------------------------:

# Application Rates for each Pesticide (source:  labels)
embrate <- mcstoc(rtriang, type="V", min=0.003, max=0.015, mode=0.015)

# Application rates with uncertainty
embfinappr <- embrate - (embrate * apprerr)

# Dermal absorption (based on registrant studies)
emabdermabs <- mcstoc(rnorm, type="U", mean=0.017776, sd=0.014668, rtrunc=TRUE, 
                      linf=0)

embinhdose <- ((nmcairdistdf * mlbrrate * embfinappr * acredis * 1000/60) + 
               (nmcairdistoc * applbrrate * embfinappr * acredis * 1000/60)) / 
    (bodyweight * 1000)

# Emamben ppe: cloth on body, gloves on hands, goggles on face
emamlderm <- (nmclowldistdf + nmcupldistdf + nmcchestdistdf + nmcbackdistdf + 
              nmcupadistdf + nmclowadistdf + 
                  nmcheaddistdf * (headnecksa/facefnecksa) - 
                  nmcheaddistdf * .1 + nmchanddistdf)
emaocderm <- (nmclowldistoc + nmcupldistoc + nmcchestdistoc + nmcbackdistoc + 
              nmcupadistoc + nmclowadistoc + 
                  (nmcheaddistoc - nmcfacedistoc * .1) + nmchanddistoc)
emaocderm <- ifelse(emaocderm>0, emaocderm, 0)
emamlderm <- ifelse(emamlderm>0, emamlderm, 0)

embdermdose <- ((emamlderm * emabdermabs * embfinappr * acredis) + 
                    (emaocderm * emabdermabs * embfinappr * acredis)) / 
    (bodyweight * 1000)
embdose <- (embdermdose + embinhdose)

efraction.exact(embdose, gam=0.95, L=0.00025 , logx=TRUE, wpnt=FALSE)
embsim <- mc(embdose)

embquant <- as.data.frame(quantile(embdose))
Percentile <- c(seq(0, 100, 1))
embquant <- t(embquant)
embquant <- as.data.frame(cbind(embquant, Percentile))
row.names(embquant) <- NULL
colnames(embquant) <- c("Median", "Mean", "LTL", "UTL", "Percentile")

#:-----------------------------------------------------------------------------:
# Dose summary 
#:-----------------------------------------------------------------------------:

summary(embsim, probs=c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), 
        lim=c(0, 0.025, 0.975, 1))

# Plots the mc object (default is median, 25th and 75th percentile, and 2.5th 
# and 97.5th percentile)
plot(embsim)

# Plots median with 2.5th and 97.5th percentiles.  
embplot <- ggplot() + 
    geom_line(data=embquant, aes(y=Percentile, x=Median)) + 
    scale_x_log10() + 
    geom_line(data=embquant, aes(y=Percentile, x=LTL), color="red") + 
    geom_line(data=embquant, aes(y=Percentile, x=UTL), color="red")
embplot

#:-----------------------------------------------------------------------------:
# Dose Response
#:-----------------------------------------------------------------------------:

# These parameters are output from the EPA BMDS software 
# (model selection was done there)

# Quantal - Linear model

# Emamectin B tremors
bckgrd2 <- mcstoc(runif, type="U", min=0, max=0)	
slope2 <- mcstoc(rnorm, type="U", mean=0.115613, sd=4.55E-05, rtrunc=TRUE, 
                 linf=-0.176931, lsup=0.408158)	

# function Emamectin B neuro (tremors - quantal)
qln2 <- bckgrd2 + (1 - bckgrd2) * (1 - exp(-slope2 * embdose))
ebneuro <- mc(qln2)
summary(ebneuro)

#:-----------------------------------------------------------------------------:
# Plots of dose - response outcome
#:-----------------------------------------------------------------------------:

plot(ebneuro)

# Extraction of  quantiles for graphic creation
embdrnquant <- as.data.frame(quantile(ebneuro))
Percentile <- c(seq(0, 100, 1))
embdrnquant <- t(embdrnquant)
embdrnquant <- as.data.frame(cbind(embdrnquant, Percentile))
row.names(embdrnquant) <- NULL
colnames(embdrnquant) <- c("Median", "Mean", "LTL", "UTL", "Percentile")

embplotneuro <- ggplot() + 
    geom_line(data=embdrnquant, aes(y=Percentile, x=Median)) + 
    scale_x_log10() + 
    geom_line(data=embdrnquant, aes(y=Percentile, x=LTL), color="red") + 
    geom_line(data=embdrnquant, aes(y=Percentile, x=UTL), color="red") + 
    xlab("Proportion of population with tremors")
embplotneuro

