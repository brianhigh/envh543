# 2-D Monte Carlo simulation of microbial exposure
Jane Pouzou  
License: CC BY-SA 4.0  

This document offers a 2-D Monte Carlo probabilistic solution in R for the 
daily microbial exposure from drinking water consumption, swimming in surface 
water and shellfish consumption for [Example 6.18](images/ex0618.png) from 
pages 215-216 of:

[Quantitative Microbial Risk Assessment, 2nd Edition](http://www.wiley.com/WileyCDA/WileyTitle/productCd-1118145291,subjectCd-CH20.html) 
by Charles N. Haas, Joan B. Rose, and Charles P. Gerba. (Wiley, 2014).


```r
# A 2-D Monte Carlo simulation for the daily microbial exposure from drinking
# water consumption, swimming in surface water, and shellfish consumption
# for Example 6.18 from Quantitative Microbial Risk Assessment, 2nd Edition
# by Charles N. Haas, Joan B. Rose, and Charles P. Gerba. (Wiley, 2014).

# Copyright (c) Jane Pouzou
# License: CC BY-SA 4.0 - https://creativecommons.org/licenses/by-sa/4.0/

# ---------------------------------------------------------------------
# Clear workspace
# ---------------------------------------------------------------------

# Clear the workspace, unless you are running in knitr context.
if (!isTRUE(getOption('knitr.in.progress'))) {
    closeAllConnections()
    rm(list = ls())
}

# ---------------------------------------------------------------------
# Set options
# ---------------------------------------------------------------------

# Set display options for use with the print() function.
options(digits = 5)

# ---------------------------------------------------------------------
# Define variables
# ---------------------------------------------------------------------

# Define the deterministic factors.
shell.viral.load <- 1
dw.viral.load <- 0.001
shell.cons <- 0.135
sw.viral.load <- 0.01
sw.frequency <- 7
seed <- 1

# Define a function to calculate risk.
Risk.fcn <- function(shell.vl, shell.cons, water.cons.L, dw.vl, sw.vl,
                     sw.daily.IR, sw.duration, sw.frequency) {
    ((shell.vl * shell.cons) + (water.cons.L * dw.vl) +
       ((sw.vl * (sw.daily.IR * sw.duration * sw.frequency)) / 365 / 1000))
}

# Define an empty matric to hold the simulation results.
Risk.mat <- matrix(as.numeric(NA), nrow = 5000, ncol = 250)

# ---------------------------------------------------------------------
# Perform a 2-D Monte Carlo simulation
# ---------------------------------------------------------------------

# Generate 250 random samples from a normal distribution to estimate 
# ingestion rate (IR) in mL of surface water. The standard deviation 
# value of 45 used here is a fictitious example.
set.seed(seed)
sw.d.IR <- rnorm(250, mean = 50, sd = 45)

# Plot the kernel density estimates for surface water.
plot(density(sw.d.IR))
```

![](ex0618prob2d_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
# Run 250 iterations of a 5000-sample simulation.
for (i in 1:250) {
    # Generate 5000 random samples from a log-normal distribution to estimate 
    # exposure from consumption of drinking water (ml/day). Divide by 1000 
    # mL/L to get consumption in liters/day.  Values for meanlog and sdlog 
    # are from the QMRA textbook (Haas et al. 2014), page 216, Table 6.30.
    set.seed(seed)
    water.cons.L <- rlnorm(n = 5000, meanlog = 7.49, sdlog = 0.407) / 1000
    
    # Plot the kernel density estimates for drinking water.
    #plot(density(water.cons.L))  # Uncomment this if you want to see these.
    
    # Sample 5000 times from a discrete distribution of swim duration with 
    # assigned probabilities of each outcome. These values are hypothetical 
    # and are not found in the text, but are defined here to provide an 
    # example of sampling from a discrete distribution.
    set.seed(seed)
    swim.duration <- sample(x = c(0.5, 1, 2, 2.6), size = 5000, 
                            replace = TRUE, prob = c(0.1, 0.1, 0.2, 0.6))
    
    # Compute 5000 daily dose simulations and store as a vector in a matrix.
    Risk.mat[,i] <- sapply(1:5000, function(j) 
        # Define a function to calculate microbial exposure risk.
        Risk.fcn(water.cons.L = water.cons.L[j],
                 sw.duration = swim.duration[j],
                 shell.vl = shell.viral.load,
                 dw.vl = dw.viral.load,
                 shell.cons = shell.cons,
                 sw.vl = sw.viral.load,
                 sw.daily.IR = sw.d.IR[i],
                 sw.frequency = sw.frequency))
}

# Plot the empirical cumulative distribution for the first iteration.
plot(ecdf(log10(Risk.mat[, 1])))

# Plot empirical cumulative distributions of additional iterations in blue.
for (j in 2:250) {
    plot(ecdf(log10(Risk.mat[, j])), col = "lightblue", add = TRUE)
}
```

![](ex0618prob2d_files/figure-html/unnamed-chunk-1-2.png)<!-- -->

```r
# ---------------------------------------------------------------------
# Repeat the simulation with mc2d
# ---------------------------------------------------------------------

# Define a function to conditionally install and load a package.
load.pkg <- function(pkg) {
    if (! suppressWarnings(require(pkg, character.only = TRUE)) ) {
        install.packages(pkg, repos = 'http://cran.r-project.org')
        library(pkg, character.only = TRUE, quietly = TRUE)
    }
}

# Load the mc2d package, installing first if necessary.
suppressMessages(load.pkg("mc2d"))
# Or use library(mc2d) and take your chances...

# Set the number of simulations in the variability dimension.
ndvar(5000)
```

```
## [1] 5000
```

```r
# Set the number of simulations in the uncertainty dimension. 
ndunc(250)
```

```
## [1] 250
```

```r
# Create mcnode objects and other variables for use with mc() and mcmodelcut().
#
# Values are from Example 6.18 from Quantitative Microbial Risk Assessment, 
# 2nd Edition by Charles N. Haas, Joan B. Rose, and Charles P. Gerba. 
# (Wiley, 2014), p. 215. Use of additional fictitious values are noted below. 
shell.vl <- mcstoc(runif, type = "V", min = 1, max = 1)
dw.vl <- mcstoc(runif, type = "V", min = 0.001, max = 0.001)
shell.cons <- mcstoc(runif, type = "V", min = 0.135, max = 0.135)
sw.vl <- mcstoc(runif, type = "V", min = 0.01, max = 0.01)
sw.frequency <- mcstoc(runif, type = "V", min = 7, max = 7)

# The standard deviation value of 45 used here is a fictitious example.
sw.daily.IR <- mcstoc(rnorm, type = "U", mean = 50, sd = 45, seed = seed, 
                      rtrunc = TRUE, linf = 0)

# Values for meanlog and sdlog are from (Haas, at al. 2014) page 216.
water.cons.L <- mcstoc(rlnorm, type = "V", meanlog = 7.49, sdlog = 0.407, 
                       seed = seed) / 1000

# The discrete distribution used here is a fictitious example.
sw.duration <- mcstoc(rempiricalD, type = "V", values = c(0.5, 1, 2, 2.6), 
                      prob = c(0.1, 0.1, 0.2, 0.6))

# There are multiple ways to run the 2D simulation depending on the 
# desired output. We will use mc() and mcmodelcut() from the mc2d package.

# Create a Monte Carlo object.
dose1 <- mc((shell.vl * shell.cons) + (water.cons.L * dw.vl) + 
                ((sw.vl * (sw.daily.IR * sw.duration * sw.frequency)
                  ) / 365 / 1000))

# Plot the Monte Carlo object.
plot(dose1)
```

![](ex0618prob2d_files/figure-html/unnamed-chunk-1-3.png)<!-- -->

```r
# Evaluate the 2-D Monte Carlo model.
o <- capture.output(dosemccut <- mcmodelcut({
    # ------------
    # First block:
    # ------------
    # Evaluate all the 0, V and U nodes.
    {
        shell.vl <- mcstoc(runif, type = "V", min = 1, max = 1)
        dw.vl <- mcstoc(runif, type = "V", min = 0.001, max = 0.001)
        shell.cons <- mcstoc(runif, type = "V", min = 0.135, max = 0.135)
        sw.vl <- mcstoc(runif, type = "V", min = 0.01, max = 0.01)
        sw.frequency <- mcstoc(runif, type = "V", min = 7, max = 7)
        
        sw.daily.IR <- mcstoc(rnorm, type = "U", mean = 50, sd = 45, 
                              seed = seed, rtrunc = TRUE, linf = 0)
        
        water.cons.L <- mcstoc(rlnorm, type = "V", meanlog = 7.49, 
                               sdlog = 0.407, seed = seed) / 1000
        
        sw.duration <- mcstoc(rempiricalD, type = "V", 
                              values = c(0.5, 1, 2, 2.6), 
                              prob = c(0.1, 0.1, 0.2, 0.6))
    }
    # -------------
    # Second block:
    # -------------
    # Evaluate all the VU nodes.
    # Leads to the mc object.
    {
        dose2 <-  ((shell.vl * shell.cons) + (water.cons.L * dw.vl) +
                     ((sw.vl * (sw.daily.IR * sw.duration * sw.frequency)
                       ) / 365 / 1000))
        
        dosemod <- mc(shell.vl, shell.cons, water.cons.L, dw.vl, sw.vl,
                      sw.daily.IR, sw.duration, sw.frequency, dose2)
    }
    # ------------
    # Third block:
    # ------------
    # Leads to a list of statistics: summary, plot, tornado
    # or any function leading to a vector (et), a list (minmax),
    # a matrix or a data.frame (summary)
    {
        list(sum = summary(dosemod), plot = plot(dosemod, draw = TRUE),
             minmax = lapply(dosemod, range))
    }
}))

# Turn off plotting if running in knitr context to save space in the report.
if (isTRUE(getOption('knitr.in.progress'))) {
    # The simulation makes a lot of plots. Don't show them.
    o <- capture.output(dev.off())
}

# Capture the text output (a progress bar) to "o" to save space in the report.
o <- capture.output(
    x <- evalmccut(dosemccut, nsv = 5000, nsu = 250, seed = seed))
print(o)
```

```
## [1] "'0' mcnode(s) built in the first block:  "                                                                     
## [2] "'V' mcnode(s) built in the first block: dw.vl shell.cons shell.vl sw.duration sw.frequency sw.vl water.cons.L "
## [3] "'U' mcnode(s) built in the first block: sw.daily.IR "                                                          
## [4] "'VU' mcnode(s) built in the first block:  "                                                                    
## [5] "The 'U' and 'VU' nodes will be sent column by column in the loop"                                              
## [6] "---------|---------|---------|---------|---------|"                                                            
## [7] "**************************************************"
```

```r
summary(x)
```

```
## shell.vl :
##       mean sd Min 2.5% 25% 50% 75% 97.5% Max  nsv Na's
## NoUnc    1  0   1    1   1   1   1     1   1 5000    0
## 
## shell.cons :
##        mean sd   Min  2.5%   25%   50%   75% 97.5%   Max  nsv Na's
## NoUnc 0.135  0 0.135 0.135 0.135 0.135 0.135 0.135 0.135 5000    0
## 
## water.cons.L :
##       mean    sd   Min  2.5%  25%  50%  75% 97.5%  Max  nsv Na's
## NoUnc 1.95 0.841 0.402 0.776 1.36 1.78 2.38  4.07 8.44 5000    0
## 
## dw.vl :
##        mean sd   Min  2.5%   25%   50%   75% 97.5%   Max  nsv Na's
## NoUnc 0.001  0 0.001 0.001 0.001 0.001 0.001 0.001 0.001 5000    0
## 
## sw.vl :
##       mean sd  Min 2.5%  25%  50%  75% 97.5%  Max  nsv Na's
## NoUnc 0.01  0 0.01 0.01 0.01 0.01 0.01  0.01 0.01 5000    0
## 
## sw.daily.IR :
##       NoVar
## 50%    57.2
## mean   61.9
## 2.5%   10.0
## 97.5% 131.3
## Nas     0.0
## 
## sw.duration :
##       mean    sd Min 2.5% 25% 50% 75% 97.5% Max  nsv Na's
## NoUnc 2.11 0.727 0.5  0.5   2 2.6 2.6   2.6 2.6 5000    0
## 
## sw.frequency :
##       mean sd Min 2.5% 25% 50% 75% 97.5% Max  nsv Na's
## NoUnc    7  0   7    7   7   7   7     7   7 5000    0
## 
## dose2 :
##        mean       sd   Min  2.5%   25%   50%   75% 97.5%   Max  nsv Na's
## 50%   0.137 0.000841 0.135 0.136 0.136 0.137 0.137 0.139 0.143 5000    0
## mean  0.137 0.000841 0.135 0.136 0.136 0.137 0.137 0.139 0.143 5000    0
## 2.5%  0.137 0.000841 0.135 0.136 0.136 0.137 0.137 0.139 0.143 5000    0
## 97.5% 0.137 0.000841 0.135 0.136 0.136 0.137 0.137 0.139 0.144 5000    0
## Nas   0.000 0.000000 0.000 0.000 0.000 0.000 0.000 0.000 0.000    0    0
```