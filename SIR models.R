
###################################################
### chunk number 1: 
###################################################
require(deSolve)                          #deSolve library needed for this computing session


###################################################
### chunk number 2: 
###################################################
sir.model <- function (t, x, params) {    #here we begin a function with three arguments
  S <- x[1]                               #create local variable S, the first element of x
  I <- x[2]                               #create local variable I
  R <- x[3]                               #create local variable R
  with(                                   #we can simplify code using "with"
    as.list(params),                   #this argument to "with" lets us use the variable names
    {                                  #the system of rate equations
      dS <- mu*(1-S)-beta*S*I
      dI <- beta*S*I-(mu+gamma)*I
      dR <- gamma*I-mu*R
      dx <- c(dS,dI,dR)                #combine results into a single vector dx
      list(dx)                         #return result as a list
    }
  )
}


###################################################
### chunk number 3: 
###################################################
times <- seq(0,10,by=1/120)              #function seq returns a sequence
params <- c(mu=1/50,beta=1000,gamma=365/13) #function c "c"ombines values into a vector
xstart <- c(S=0.06,I=0.001,R=0.939)      #initial conditions


###################################################
### chunk number 4: 
###################################################
out <- as.data.frame(lsoda(xstart,times,sir.model,params))  #result stored in dataframe


###################################################
### chunk number 5: 
###################################################
op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))                  #set graphical parameters
plot(I~time,data=out,type='l',log='y')                      #plot the I variable as a line
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)                  #re-set graphical parameters
plot(I~S,data=out,type='p',log='xy',yaxt='n',xlab='S',cex=0.5)  #plot phase portrait
par(op)                                                     #re-set graphical parameters


###################################################
### chunk number 6: 
###################################################
seasonal.sir.model <- function (t, x, params) {  #function to return time-dependent rates
  with(
    as.list(c(x,params)),
    {
      beta <- beta0*(1+beta1*cos(2*pi*t))     #first determine time-dependent beta
      dS <- mu*(1-S)-beta*S*I                 #the system of rate equations
      dI <- beta*S*I-(mu+gamma)*I
      dR <- gamma*I-mu*R
      dx <- c(dS,dI,dR)                       #store result
      list(dx)                                #and return as a list
    }
  )
}


###################################################
### chunk number 7: 
###################################################
times <- seq(0,100,by=1/120)                     #times at which to obtain solution
params <- c(mu=1/50,beta0=1000,beta1=0.4,gamma=365/13)   #parameters
xstart <- c(S=0.06,I=0.001,R=0.939)              #initial conditions
out <- as.data.frame(lsoda(xstart,times,seasonal.sir.model,params,rtol=1e-12,hmax=1/120)) #solve
op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))       #set graphical parameters
plot(I~time,data=out,type='l',log='y',subset=time>=90) #plot
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)       #reset graphical parameters
plot(I~S,data=out,type='p',log='xy',subset=time>=50,yaxt='n',xlab='S',cex=0.5) #plot
text(0.02,1e-02,cex=0.7, "last 50 yr of simulation")  #annotate graph
par(op)                                          #reset graphical parameters


###################################################
### chunk number 8: 
###################################################
seir.model <- function (t, x, params) {    #here we begin a function with three arguments
  S <- x[1]                               #create local variable S, the first element of x
  E <- x[2]                               #create local variable E
  I <- x[3]                               #create local variable I
  R <- x[4]                               #create local variable R
  with(                                   #we can simplify code using "with"
    as.list(params),                   #this argument to "with" lets us use the variable names
    {                                  #the system of rate equations
      dS <- mu*(1-S)-beta*S*I
      dE <- beta*S*I-(sigma+mu)*E
      dI <- sigma*E-(mu+gamma)*I
      dR <- gamma*I-mu*R
      dx <- c(dS,dE,dI,dR)             #combine results into a single vector dx
      list(dx)                         #return result as a list
    }
  )
}


###################################################
### chunk number 9: 
###################################################
times <- seq(0,10,by=1/120)              #function seq returns a sequence
params <- c(mu=1/50,sigma=50,beta=50,gamma=365/104) #function c "c"ombines values into a vector
xstart <- c(S=0.45,E=0.10,I=0.050,R=0.400)      #initial conditions


###################################################
### chunk number 10: 
###################################################
out <- as.data.frame(lsoda(xstart,times,seir.model,params))  #result stored in dataframe


###################################################
### chunk number 11: 
###################################################
op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))                  #set graphical parameters
plot(I~time,data=out,type='l',log='y')                      #plot the I variable as a line
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)                  #re-set graphical parameters
plot(I~S,data=out,type='p',log='xy',yaxt='n',xlab='S',cex=0.5)  #plot phase portrait
par(op)                                                     #re-set graphical parameters
