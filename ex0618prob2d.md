# 2-D Monte Carlo simulation of microbial exposure
Jane Pouzou and Brian High  
![CC BY-SA 4.0](cc_by-sa_4.png)  

## Introduction

This document offers three 2-D Monte Carlo probabilistic solutions in R for the 
daily microbial exposure from drinking water consumption, swimming in surface 
water and shellfish consumption for [Example 6.18](images/ex0618.png) from 
pages 215-216 of:

[Quantitative Microbial Risk Assessment, 2nd Edition](http://www.wiley.com/WileyCDA/WileyTitle/productCd-1118145291,subjectCd-CH20.html) 
by Charles N. Haas, Joan B. Rose, and Charles P. Gerba. (Wiley, 2014).

The first solution is coded manually with functionality built into "base" R.

The second and third solutions use the _mc2d_ package.

## Set global options


```r
# Set display options for use with the print() function.
options(digits = 5)

# Set knitr options for use when rendering this document.
library("knitr")
opts_chunk$set(cache=TRUE, message=FALSE)
```

## Perform a 2-D Monte Carlo simulation

### Define variables

1. Define the deterministic factors from the microbial exposure example. 
2. Define the _seed_ used for random sampling to enable reproducible results. 
   This variable will be used later to set the seed using the `set.seed()`
   function before random sampling. 
3. Define the number of random samples to be drawn from probability 
   distrubutions in the two dimensions (variability and uncertainty) of the 
   simulation.


```r
# Values from Example 6.18 from Quantitative Microbial Risk Assessment, 
# 2nd Edition by Charles N. Haas, Joan B. Rose, and Charles P. Gerba. 
# (Wiley, 2014), pp. 215-216.
shellfish.viral.load  <- 1        # Shellfish viral loading (viruses/g)
dw.viral.load         <- 0.001    # Drinking water viral loading (viruses/L)
shellfish.cons.g      <- 0.135    # Shellfish consumption (viruses/day)
sw.viral.load         <- 0.1      # Swimming in viral loading (viruses/L)
sw.frequency          <- 7        # Swimming frequency (swims/year)

# Define an integer to use when setting the seed of the random number generator.
seed <- 1

# Define the number of samples for the variability and uncertainty dimensions.
nsv <- 5000
nsu <- 250
```

### Sample from distributions reflecting uncertainty 

In our previous [1-D Monte Carlo simulation](ex0618prob.md), we used a point 
estimate for the ingestion rate of surface water while swimming. Now, for our 
2-D simulation, we will represent the uncertainty around this estimate with a 
normal distribution.

For the uncertainty dimension, generate 250 random samples from a normal 
distribution to estimate ingestion rate (IR) in mL of surface water while
swimming. 


```r
# The point estimate of 50 mL/hr for ingestion of surface water while swimming 
# is taken from the QMRA textbook (Haas et al. 2014), page 215. We set it here
# as the mean. The standard deviation value of 45 used here is fictitious.
set.seed(seed)
sw.d.IR <- rnorm(nsu, mean = 50, sd = 45)
```

Plot the kernel density estimates for surface water ingestion rate.


```r
plot(density(sw.d.IR))
```

![](ex0618prob2d_files/figure-html/kernel-density-plot-1.png)<!-- -->

### Define exposure risk function

Define a risk function to calculate estimated microbial exposure.


```r
Risk.fcn <- function(shellfish.vl, shellfish.cons.g, water.cons.L, dw.vl, 
                     sw.vl, sw.daily.IR, sw.duration, sw.frequency) {
    ((shellfish.vl * shellfish.cons.g) + (water.cons.L * dw.vl) + 
         ((sw.vl * (sw.daily.IR * sw.duration * sw.frequency)) / 365 / 1000))
}
```

### Run a 2-D MC simulation

Run the simulation using a nested loop structure of 5000 iterations (in the 
variability dimension) performed by the `sapply()` function for each of 250 
iterations (in the uncertainty dimension) of a `for()` loop.

Specifically, within the `for()` loop:

- Sample from the distributions reflecting variability:
    1. Drinking water consumption in liters/day (log-normal distribution)
    2. Swimming duration in hours (discrete distribution)
- Compute 5000 daily exposure estimates and store as a vector in a matrix.

The results of each iteration will be accumulated in a 5000 row by 250 column 
matrix. 


```r
# Define an empty matrix to hold the simulation results.
Risk.mat <- matrix(as.numeric(NA), nrow = nsv, ncol = nsu)

# Loop through the uncertainty dimension.
for (j in 1:nsu) {
    # 1. Generate 5000 random samples from a log-normal distribution to estimate 
    #    exposure from consumption of drinking water (ml/day). Divide by 1000 
    #    mL/L to get consumption in liters/day.  Values for meanlog and sdlog 
    #    are from the QMRA textbook (Haas et al. 2014), page 216, Table 6.30.
    set.seed(seed)
    water.cons.L <- rlnorm(n = nsv, meanlog = 7.49, sdlog = 0.407) / 1000
    
    # 2. Sample 5000 times from a discrete distribution of swim duration with 
    #    assigned probabilities of each outcome. These values are hypothetical 
    #    and are not found in the text, but are defined here to provide an 
    #    example of sampling from a discrete distribution.
    set.seed(seed)
    swim.duration <- sample(x = c(0.5, 1, 2, 2.6), size = nsv, 
                            replace = TRUE, prob = c(0.1, 0.1, 0.2, 0.6))
    
    # Loop through the variability dimension.
    # Compute 5000 daily simulations and store as a vector in a matrix.
    Risk.mat[, j] <- sapply(1:nsv, function(i) 
        # Define a function to calculate esimated microbial exposure risk.
        Risk.fcn(water.cons.L = water.cons.L[i],
                 sw.duration = swim.duration[i],
                 shellfish.vl = shellfish.viral.load,
                 dw.vl = dw.viral.load,
                 shellfish.cons.g = shellfish.cons.g,
                 sw.vl = sw.viral.load,
                 sw.daily.IR = sw.d.IR[j],
                 sw.frequency = sw.frequency))
}
```

### Summarize results

Since this is a 2-D simulation, we can report the uncertainty in our point 
estimate of the microbial exposure risk.

We will use the 50th percentile (median) of the means of each iteration in the
uncertainty dimension as our point estimate. The 2.5th and 97.5th percentiles 
are used to establish a 95% confidence interval (CI95).


```r
# Report the median of the means with a 95% confidence interval (CI95).
mean.risk <- sapply(1:nsu, function(j) mean(Risk.mat[, j]))
quantile(mean.risk, probs = seq(0, 1, 0.025))[c("2.5%", "50%", "97.5%")]
```

```
##    2.5%     50%   97.5% 
## 0.13683 0.13714 0.13752
```

Next, we build up an empirical cumulative distribution plot of all of the 
loop iterations, by overlaying each new curve upon a plot of the first 
iteration.


```r
# Plot the empirical cumulative distribution for the first iteration.
plot(ecdf(Risk.mat[, 1]), col = "#ADD8E605", main = 'ecdf(Risk.mat)')

# Plot empirical cumulative distributions for additional iterations in blue.
for (j in 2:nsu) {
    plot(ecdf(Risk.mat[, j]), col = "#ADD8E605", add = TRUE)
}
```

![](ex0618prob2d_files/figure-html/ecdf-plot-risk-mat-1.png)<!-- -->

## Repeat the simulation with mc2d

We will demonstrate two alternative ways to run 2-dimensional Monte Carlo 
simulation using the 
[mc2d](https://cran.r-project.org/web/packages/mc2d/index.html) package.

For our first alternative, we will use 
[mcmodel](http://www.inside-r.org/packages/cran/mc2d/docs/mcmodel) and 
[evalmcmod](http://www.inside-r.org/packages/cran/mc2d/docs/evalmcmod) from 
the [mc2d](https://cran.r-project.org/web/packages/mc2d/index.html) package. 
This is a simple way to implement the 2-D simulation because we just need to
define the model with `mcmodel()` and evaluate it with `evalmcmod()`.

### Load packages

Make sure the _mc2d_ package is installed before trying to use it. If we just 
tried to load it with `library()` or `require()` before it had been installed, 
our program would end with errors. (That would not be very friendly.)


```r
# Define a function to conditionally install and load a package.
load.pkg <- function(pkg) {
    if (! suppressWarnings(require(pkg, character.only = TRUE)) ) {
        install.packages(pkg, repos = 'http://cran.r-project.org')
        library(pkg, character.only = TRUE, quietly = TRUE)
    }
}

# Load required packages, installing first if necessary.
suppressMessages(load.pkg("mc2d"))   # Or just use: library(mc2d)
```

### Define variables

Set the number of simulations for the variability and uncertainty dimensions. 


```r
ndvar(5000)  # Variability
```

```
## [1] 5000
```

```r
ndunc(250)   # Uncertainty
```

```
## [1] 250
```

The _mc2d_ functions will set `nsv` to `ndvar()` and `nsu` to `ndunc()` by 
default. Alternatively, we could supply these values to the _mc2d_ model
functions when we call them.

Define the variables we will use to set the `seed` for random sampling and the 
number of `digits` for `print()` statements.


```r
seed <- 1
digits <- 5
```

We will use this variable to explicitly set the seed with the various
_mc2d_ functions. Another approach would be to only set it through the
model evaluation functions, or not at all. Since we want to do our best
to provide the most reproducible results, we set the seed explicitly.

### Define exposure model

Within the [mcmodel](http://www.inside-r.org/packages/cran/mc2d/docs/mcmodel)
function, use [mcstoc](http://www.inside-r.org/packages/cran/mc2d/docs/mcstoc)
to define "mc nodes" for each component of the model. 

For each stochastic variable, use the `mcstoc()` function to create a 
"mc node", supply a probability function, the node type as "V" for variability 
or "U" for uncertainty, and any additional parameters to be passed to the 
probability function. 

For deterministic variables, create nodes with the `mcdata()` function 
using the "0" type.

We will model the deterministic factors as uniform probablity distributions.

The last statement makes an "mc" object using the 
[mc](http://www.inside-r.org/packages/cran/mc2d/docs/mc) function.


```r
# Define an exposure model for evaluation by evalmcmod().
expo.mod1 <- mcmodel({
    # Values from Example 6.18 from Quantitative Microbial Risk Assessment, 
    # 2nd Edition by Charles N. Haas, Joan B. Rose, and Charles P. Gerba. 
    # (Wiley, 2014), pp. 215-216. Other fictitious values are noted below.
        
    # Shellfish viral loading (viruses/g):
    shellfish.vl <- mcdata(1, type = "0")
    
    # Shellfish consumption (g/day):
    shellfish.cons.g <- mcdata(0.135, type = "0")
    
    # Drinking water viral loading (viruses/L):
    dw.vl <- mcdata(0.001, type = "0")
    
    # Drinking water consumption (L/day):
    dw.cons.L <- mcstoc(rlnorm, type = "V", seed = seed, 
                        meanlog = 7.49, sdlog = 0.407) / 1000
    
    # Swimming in surface water viral loading (viruses/L):
    sw.vl <- mcdata(0.1, type = "0")
    
    # Swimming daily ingestion rate (mL/hour): fictitious sd = 45
    sw.daily.IR <- mcstoc(rnorm, type = "U", seed = seed, 
                          mean = 50, sd = 45, rtrunc = TRUE, linf = 0)
    
    # Swimming duration (hours): fictitious discrete distribution
    sw.duration <- mcstoc(rempiricalD, type = "V", seed = seed, 
                          values = c(0.5, 1, 2, 2.6), 
                          prob = c(0.1, 0.1, 0.2, 0.6))
    
    # Swimming frequency (swims/year):
    sw.frequency <- mcdata(7, type = "0")
    
    # Estimate the exposure using the 0, V and U nodes to create a VU node.
    expo.mc1 <- (shellfish.vl * shellfish.cons.g) + (dw.vl * dw.cons.L) + 
           ((sw.vl * (sw.daily.IR * sw.duration * sw.frequency)) / 365 / 1000)
    
    # Build a mc model from all of the mcnode objects.
    mc(shellfish.vl, shellfish.cons.g, dw.vl, dw.cons.L, sw.vl, sw.daily.IR, 
       sw.duration, sw.frequency, expo.mc1) 
})
```

### Evaluate the model

Evaluate the model with 5000 iterations in the variability dimension and 250 
iterations in the uncertainty dimension, as set previously.


```r
expo.ev1 <- evalmcmod(expo.mod1, seed = seed)
print(expo.ev1, digits = digits)
```

```
##   node    mode  nsv nsu nva variate     min   mean  median     max Nas
## 1      numeric 5000 250   1       1 0.13541 0.1372 0.13704 0.14425   0
##   type outm
## 1   VU each
```

### Summarize results

Print a summary and a plot of the evaluation results (`expo.ev1`). Report the 
median of the means with a 95% confidence interval (CI95). Plot the empirical 
cumulative distribution function (ecdf) of the exposure model (`expo.mc1`) 
estimates.


```r
# Summarize the results.
summary(expo.ev1)
```

```
##  :
##         mean       sd   Min  2.5%   25%   50%   75% 97.5%   Max  nsv Na's
## median 0.137 0.000845 0.136 0.136 0.137 0.137 0.138 0.139 0.144 5000    0
## mean   0.137 0.000847 0.136 0.136 0.137 0.137 0.138 0.139 0.144 5000    0
## 2.5%   0.137 0.000841 0.135 0.136 0.136 0.137 0.137 0.139 0.143 5000    0
## 97.5%  0.137 0.000861 0.136 0.136 0.137 0.137 0.138 0.140 0.144 5000    0
```

```r
# Plot the results.
plot(expo.ev1)
```

![](ex0618prob2d_files/figure-html/results-ev1-1.png)<!-- -->

```r
# Report the median of the means with a 95% confidence interval (CI95).
#mean.risk1 <- sapply(1:ndunc(), function(j) mean(expo.ev1$expo.mc1[, j, ]))
#quantile(mean.risk1, probs = seq(0, 1, 0.025))[c("2.5%", "50%", "97.5%")]

# Generate an "ecdf" plot.
#plot(expo.ev1$expo.mc1)
```

## Repeat 2-D simulation again with a loop

This time we will use another pair of functions from the 
[mc2d](https://cran.r-project.org/web/packages/mc2d/index.html) 
package, `mcmodelcut()` and `evalmccut()`, to get a different 
style of summary output. They implement the uncertainty dimension of the 
2-D simulation with a processing _loop_.

By using these [mccut](http://www.inside-r.org/packages/cran/mc2d/docs/mccut) 
functions, we will also be able to conserve system memory, at the cost of 
taking more time to perform the simulation. This approach allows the evaluation 
of very high dimensional models on relatively modest computer systems.

### Load packages


```r
# Load packages.
library(mc2d)
```

### Define variables

Define the variables we will use to set the `seed` for random sampling and the 
number of `digits` for `print()` statements.


```r
seed <- 1
digits <- 5
```

### Define exposure model

The `mcmodelcut()` function expects its input in three blocks:

```
mcmodelcut({ 
    # Block 1: Evaluate all of the 0, V and U mc node objects.
    { 
        # This block is evaluated once before the first loop (step 1).
    } 
    # Block 2: Evaluate all of the VU nodes. Last statement makes an mc object.
    { 
        # This block is evaluated using nsu = 1 (step 2).
    } 
    # Block 3: Build a list of statistics refering to the mc object.
    { 
        # This block is repeated nsu times (step 3).
    } 
})
```

... where `nsu` is the number of simulations for uncertainty used in the 
evaluation.


```r
# Build a mcmodelcut object in three blocks for evaluation by evalmccut().
expo.mcmcut <- mcmodelcut({
    # Block 1: Evaluate all of the 0, V and U mc node objects.
    {
        # Values from Example 6.18 from Quantitative Microbial Risk Assessment, 
        # 2nd Edition by Charles N. Haas, Joan B. Rose, and Charles P. Gerba. 
        # (Wiley, 2014), pp. 215-216. Other fictitious values are noted below.
        
        # Shellfish viral loading (viruses/g):
        shellfish.vl <- mcdata(1, type = "0")
        
        # Shellfish consumption (g/day):
        shellfish.cons.g <- mcdata(0.135, type = "0")
        
        # Drinking water viral loading (viruses/L):
        dw.vl <- mcdata(0.001, type = "0")
        
        # Drinking water consumption (L/day):
        dw.cons.L <- mcstoc(rlnorm, type = "V", seed = seed, 
                            meanlog = 7.49, sdlog = 0.407) / 1000
        
        # Swimming in surface water viral loading (viruses/L):
        sw.vl <- mcdata(0.1, type = "0")
        
        # Swimming daily ingestion rate (mL/hour): fictitious sd = 45
        sw.daily.IR <- mcstoc(rnorm, type = "U", seed = seed, 
                              mean = 50, sd = 45, rtrunc = TRUE, linf = 0)
        
        # Swimming duration (hours): fictitious discrete distribution
        sw.duration <- mcstoc(rempiricalD, type = "V", seed = seed, 
                              values = c(0.5, 1, 2, 2.6), 
                              prob = c(0.1, 0.1, 0.2, 0.6))
        
        # Swimming frequency (swims/year):
        sw.frequency <- mcdata(7, type = "0")
    }
    
    # Block2: Evaluate all of the VU nodes. Last statement makes an mc object.
    {
        # Estimate the exposure using the 0, V and U nodes to create a VU node.
        expo.mc2 <- (shellfish.vl * shellfish.cons.g) + (dw.vl * dw.cons.L) + 
            ((sw.vl * (sw.daily.IR * sw.duration * sw.frequency)) / 365 / 1000)
   
        # Build a mc model from all of the mcnode objects.
        expo.mod2 <- mc(shellfish.vl, shellfish.cons.g, dw.vl, dw.cons.L, 
                        sw.vl, sw.daily.IR, sw.duration, sw.frequency, expo.mc2)
    }

    # Block 3: Build a list of statistics refering to the mc object.
    {
        list(sum = summary(expo.mod2), plot = plot(expo.mod2, draw = FALSE))
    }

})
```

```
## The following expression will be evaluated only once :
## {
##     shellfish.vl <- mcdata(1, type = "0")
##     shellfish.cons.g <- mcdata(0.135, type = "0")
##     dw.vl <- mcdata(0.001, type = "0")
##     dw.cons.L <- mcstoc(rlnorm, type = "V", seed = seed, meanlog = 7.49, 
##         sdlog = 0.407)/1000
##     sw.vl <- mcdata(0.1, type = "0")
##     sw.daily.IR <- mcstoc(rnorm, type = "U", seed = seed, mean = 50, 
##         sd = 45, rtrunc = TRUE, linf = 0)
##     sw.duration <- mcstoc(rempiricalD, type = "V", seed = seed, 
##         values = c(0.5, 1, 2, 2.6), prob = c(0.1, 0.1, 0.2, 0.6))
##     sw.frequency <- mcdata(7, type = "0")
## }
## The mc object is named:  expo.mod2
```

### Evaluate the model

Evaluate the model with 5000 iterations in the variability dimension and 
250 iterations in the uncertainty dimesion. Save the evaluation results as 
`expo.ev2`. 

Since the `evalmccut()` function produces a lot of text output  
that we do not want in our report, we capture the text output with the 
`capture.output()` function and print a summary when finished.


```r
capture.output(expo.ev2 <- evalmccut(expo.mcmcut, seed = seed))
```

```
## [1] "'0' mcnode(s) built in the first block: dw.vl shellfish.cons.g shellfish.vl sw.frequency sw.vl "
## [2] "'V' mcnode(s) built in the first block: dw.cons.L sw.duration "                                 
## [3] "'U' mcnode(s) built in the first block: sw.daily.IR "                                           
## [4] "'VU' mcnode(s) built in the first block:  "                                                     
## [5] "The 'U' and 'VU' nodes will be sent column by column in the loop"                               
## [6] "---------|---------|---------|---------|---------|"                                             
## [7] "**************************************************"
```

### Summarize results

Print the accumulated statistics with `summary()` and `plot()`. Report the 
median of the means with a 95% confidence interval (CI95).


```r
# Print a summary and a plot for the mccut object.
summary(expo.ev2)
```

```
## shellfish.vl :
##       NoVar
## NoInc     1
## 
## shellfish.cons.g :
##       NoVar
## NoInc 0.135
## 
## dw.vl :
##       NoVar
## NoInc 0.001
## 
## dw.cons.L :
##       mean    sd   Min  2.5%  25%  50%  75% 97.5%  Max  nsv Na's
## NoUnc 1.95 0.841 0.402 0.776 1.36 1.78 2.38  4.07 8.44 5000    0
## 
## sw.vl :
##       NoVar
## NoInc   0.1
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
## NoUnc  2.1 0.744 0.5  0.5   2 2.6 2.6   2.6 2.6 5000    0
## 
## sw.frequency :
##       NoVar
## NoInc     7
## 
## expo.mc2 :
##        mean       sd   Min  2.5%   25%   50%   75% 97.5%   Max  nsv Na's
## 50%   0.137 0.000845 0.136 0.136 0.137 0.137 0.138 0.139 0.144 5000    0
## mean  0.137 0.000847 0.136 0.136 0.137 0.137 0.138 0.139 0.144 5000    0
## 2.5%  0.137 0.000841 0.135 0.136 0.136 0.137 0.137 0.139 0.143 5000    0
## 97.5% 0.137 0.000861 0.136 0.136 0.137 0.137 0.138 0.140 0.144 5000    0
## Nas   0.000 0.000000 0.000 0.000 0.000 0.000 0.000 0.000 0.000    0    0
```

```r
plot(expo.ev2)
```

![](ex0618prob2d_files/figure-html/results-ev2-1.png)<!-- -->

```r
# Report the median of the means with a 95% confidence interval (CI95).
mean.risk2 <- expo.ev2$sum$expo.mc2[, , "mean"]
quantile(mean.risk2, probs = seq(0, 1, 0.025))[c("2.5%", "50%", "97.5%")]
```

```
##    2.5%     50%   97.5% 
## 0.13699 0.13718 0.13748
```

Plot the empirical cumulative distribution function (ecdf) for the estimated 
exposure. To match the format of the `ecdf` plot made earlier, we will 
construct the plot elements manually.


```r
expo.x <- expo.ev2$plot$expo.mc2
expo.l <- length(expo.x)
expo.y <- 1:expo.l/expo.l
plot(expo.x, expo.y, pch = 20, cex = 0.1, col = '#ADD8E605', 
     main='ecdf(expo.mc2)', ylab = 'Fn(x)', xlab = 'x')
abline(h = 0, col = "gray", lty = 2, lwd = 2)
abline(h = 1, col = "gray", lty = 2, lwd = 2)
```

![](ex0618prob2d_files/figure-html/results-ev2-ecdf-plot-1.png)<!-- -->
