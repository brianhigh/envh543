

# A probabilistic solution in R for the daily microbial exposure from drinking
# water consumption for Example 6.18 from Quantitative Microbial Risk 
# Assessment, 2nd Edition by Charles N. Haas, Joan B. Rose, and 
# Charles P. Gerba. (Wiley, 2014).

# Copyright (c) Jane Pouzou
# License: CC BY-SA 4.0 - https://creativecommons.org/licenses/by-sa/4.0/


# Set display options for use with the print() function.
options(digits = 5)

# Here are the deterministic factors.
shell.viral.load <- 1
dw.viral.load <- 0.001
sw.viral.load <- 0.01
sw.frequency <- 7


# Define a function to calculate risk.
Risk.fcn <- function(shell.vl, shell.cons, water.cons.L, dw.vl, sw.vl, 
                     sw.daily.IR, sw.duration, sw.frequency) {
  ((shell.vl * shell.cons) + (water.cons.L * dw.vl) + 
     ((sw.vl * (sw.daily.IR * sw.duration * sw.frequency)) / 365 / 1000))
}


Risk.mat <- matrix(as.numeric(NA), nrow = 5000, ncol = 250)
set.seed(1)
sw.d.IR <-rnorm(250, mean=50, sd=45)
shell.con  <- rlnorm(n = 250, meanlog = -2.0025, sdlog = 0.02)

for(i in 1:250) {
  
  set.seed(1)
  water.cons.L <- rlnorm(n = 5000, meanlog = 7.49, sdlog = 0.407) / 1000
  set.seed(1)
  swim.duration <- sample(x = c(0.5, 1, 2, 2.6), size = 5000, replace = TRUE, prob = c(0.1, 0.1, 0.2, 0.6))
  sw.daily.IR<-sw.d.IR[i]*1
  shell.cons<-shell.con[i]*1
  
  Risk.mat[ ,i] <-sapply(1:5000, function(j) Risk.fcn(water.cons.L = water.cons.L[j], 
                                                      sw.duration = swim.duration[j], 
                                                      shell.vl = shell.viral.load, 
                                                      dw.vl = dw.viral.load, 
                                                      shell.cons = shell.cons, 
                                                      sw.vl = sw.viral.load, 
                                                      sw.daily.IR = sw.daily.IR, 
                                                      sw.frequency = sw.frequency))
  
  
  
}



# Plot the results of the simulation.
plot(ecdf(Risk.mat[, 1]), col="red", add=TRUE)
for(j in 2:250) {
  plot(ecdf((Risk.mat[, j])), col="lightblue", add=TRUE, xlim=c(0, 0.17))
}


####now with mc2d
require(mc2d)

ndvar(5000)
ndunc(250)

shell.vl <- mcstoc(runif, type="V", min=1, max=1)
dw.vl <- mcstoc(runif, type="V", min=0.001, max=0.001)
shell.cons <- mcstoc(rlnorm, type="U", meanlog = -2.0025, sdlog = 0.02, seed=1)
sw.vl <- mcstoc(runif, type="V", min=0.01, max=0.01)
sw.frequency <- mcstoc(runif, type="V", min=7, max=7)
sw.daily.IR <-mcstoc(rnorm, type="U", mean=50, sd=45, seed=1, rtrunc=TRUE, linf=0)
water.cons.L<-mcstoc(rlnorm, type="V", meanlog = 7.49, sdlog = 0.407, seed=1)/ 1000
sw.duration<-mcstoc(rempiricalD, type="V", values=c(0.5, 1, 2, 2.6), prob=c(0.1, 0.1, 0.2, 0.6))

######There are multiple ways to run the 2D simulation depending on the desired output.
dose1<-mc((shell.vl * shell.cons) + (water.cons.L * dw.vl) + 
            ((sw.vl * (sw.daily.IR * sw.duration * sw.frequency)) / 365 / 1000))

plot(dose1)

dosemccut<-mcmodelcut({
## First block:
## Evaluates all the 0, V and U nodes.
{ shell.vl <- mcstoc(runif, type="V", min=1, max=1)
  dw.vl <- mcstoc(runif, type="V", min=0.001, max=0.001)
  shell.cons <- mcstoc(runif, type="V", min=0.135, max=0.135)
  sw.vl <- mcstoc(runif, type="V", min=0.01, max=0.01)
  sw.frequency <- mcstoc(runif, type="V", min=7, max=7)
  
  sw.daily.IR <-mcstoc(rnorm, type="U", mean=50, sd=45, seed=1, rtrunc=TRUE, linf=0)
  
  water.cons.L<-mcstoc(rlnorm, type="V", meanlog = 7.49, sdlog = 0.407, seed=1)/ 1000
  
  sw.duration<-mcstoc(rempiricalD, type="V", values=c(0.5, 1, 2, 2.6), prob=c(0.1, 0.1, 0.2, 0.6))
  
}
## Second block:
## Evaluates all the VU nodes
## Leads to the mc object. 
{
 
  dose2 <-  ((shell.vl * shell.cons) + (water.cons.L * dw.vl) + 
              ((sw.vl * (sw.daily.IR * sw.duration * sw.frequency)) / 365 / 1000))
 
  dosemod <- mc(shell.vl, shell.cons, water.cons.L, dw.vl, sw.vl,
                sw.daily.IR, sw.duration, sw.frequency, dose2)
}
## Third block:
## Leads to a list of statistics: summary, plot, tornado
## or any function leading to a vector (et), a list (minmax), 
## a matrix or a data.frame (summary)
{
  list(
    sum = summary(dosemod), 
    plot = plot(dosemod, draw=TRUE), 
    minmax = lapply(dosemod, range)
  )
}
})

x <- evalmccut(dosemccut, nsv = 5000, nsu = 250, seed = 1)
summary(x)

