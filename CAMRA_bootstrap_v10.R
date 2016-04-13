#==========================================================================|
#-----------------------------------HEADER---------------------------------|
#==========================================================================|
#CAMRA_bootstrap_v6.r is an R statistical programming package		   |
#source code for initializing the bootstrap technique in R for		   |
#generation of confidence intervals for both the exponential and           |
#beta Poisson dose response models for risk estimation of                  |
#exposure to pathogens                                                     |
#Written for CAMRA by: Mark H. Weir, Timothy A. Bartrand and               |
#Charles N. Haas, Department of the Civil Architectural and Environmental  |
#Engineering, Drexel University                                            |
#==========================================================================|

#==========================================================================|
#--Draw in dose response data and assign the needed values from the data---|
#==========================================================================|
DR_Data <- read.table("RoV_dr_all.txt",header=TRUE)				  
	dose <- DR_Data$dose									  
	positive <- DR_Data$positive_response									  
	negative <- DR_Data$negative_response									   
	oprob <- positive/(positive+negative)									   
#==========================================================================|

#==========================================================================|
#-----------Load required libraries to perform MLE and bootstrap-----------|
#==========================================================================|
	require(stats4)									         
	require(boot)	
	require(car)								        
#==========================================================================|

dose = DR_Data$dose
positive = DR_Data$positive_response
negative = DR_Data$negative_response
ni = positive+negative

xi = log(dose)
xbar = sum(ni*xi)/sum(ni)
pbar = sum(positive)/sum(ni)

Zca = sum((xi-xbar)*positive)/sqrt((pbar*(1-pbar))*sum(ni*(xi-xbar)^2))
if (Zca>1.644)
{

#==========================================================================|
#----------------------Define the dose response models.--------------------|
#----expntl.dr ---> exponential model and bp.dr ---> beta Poisson model----|
#==========================================================================|
	expntl.dr <- function(k,dose) 1 - exp(-k*dose)    
	bp.dr <- function(alpha,N50,dose) 1-(1+(dose/N50)*(2^(1/alpha)-1))^(-alpha)
#==========================================================================|

#==========================================================================|
#---------Define functions for deviances of dose-response model------------| 
#---------------------fits to experimental data----------------------------|
#----deviance.expntl --> deviance of the exponential model to the data ----|
#--------dev.bp --> deviance of the beta Poisson model to the data---------|
#================================ Line 44 =================================|											        
deviance.expntl <- function(obspos, obsneg, logk, dose)                   
	{                                                                         
      eps = 1e-15;       #ensures that the function will not divide by zero
      k = exp(logk)
      obsf = obspos/(obspos + obsneg);
      pred = expntl.dr(k,dose);
      y1 = sum(obspos*log(pred/(obsf+eps)));
      y2 = sum(obsneg*log((1-pred+eps)/(1-obsf+eps)))
      return(-1*(y1+y2))
	}

dev.bp <- function (obspos, obsneg, logalpha, logN50, dose)
	{
	eps = 1e-15;
	alpha = exp(logalpha)
	N50 = exp(logN50)
	obsf = obspos/(obspos+obsneg);
	pred = bp.dr(alpha,N50,dose);
	y1 = sum(obspos*log(pred/(obsf+eps)));
        y2 = sum(obsneg*log((1-pred+eps)/(1-obsf+eps)))
        return(-1*(y1+y2))
	}
#=============================== Line 67 ==================================|

#==========================================================================|
#------First run an MLE routine to obtain intial fitting estimates---------|
#------------for both exponential and beta Poisson models------------------|
#==========================================================================|
      
results2<-mle(deviance.expntl, start=list(logk=-11), 
    method = 'BFGS',fixed = list(obspos=positive,obsneg=negative,dose=dose))

	EXP_MLE <- matrix(ncol=3, nrow=1)
	YEget <- logLik(results2)
	YE <- -2*YEget[[1]]
	j<-coef(results2)
	logk=j["logk"]
	k = exp(logk)
	Exp_ID50=(log(-1/(0.5-1))/k)
	
	EXP_MLE[,1] = YE; EXP_MLE[,2] = k; EXP_MLE[,3]=Exp_ID50;

	colnames(EXP_MLE)<-c("Minimized Deviance","k Parameter","50% Probability Dose (Lethal or Infectious)")
	rownames(EXP_MLE)<-("")
	
	write.csv(EXP_MLE, file="Exponential_MLE_output.csv")

print("==========================================================")
print(" MLE results for the exponential dose response model")
	print(EXP_MLE)
print("==========================================================")


resultsbp <- mle(dev.bp,start=list(logalpha=-0.2, logN50=10),
	method='BFGS',fixed=list(obspos=positive,obsneg=negative,dose=dose))

	BP_MLE <- matrix(ncol=4,nrow=1)
	YBPget <- logLik(resultsbp)
	YBP <- -2*YBPget[[1]]
	jb<-coef(resultsbp)
	bpalpha = exp(jb["logalpha"]); bpN50 = exp(jb["logN50"])
	logalpha=jb["logalpha"]
	logN50=jb["logN50"]
	BP_ID50=abs(((0.5^(-1/bpalpha)-1)*bpN50)/(2^(1/bpalpha)-1))	

	BP_MLE[,1] <- YBP;
	BP_MLE[,2] <- bpalpha;
	BP_MLE[,3] <- bpN50;
	BP_MLE[,4] <- BP_ID50;

	colnames(BP_MLE)<- c("Minimized Deviance","alpha","N50","LD50 or ID50"); 
	rownames(BP_MLE)<-(" ")

	write.csv(BP_MLE, file="betPoisson_MLE_output.csv")

print("===========================================================")
print(" MLE results for the beta Poisson dose response model")
	print(BP_MLE)
print("===========================================================")

#============================================Line 125===========================================|

print("======================================================================================")
print("   Evaluate the Goodness of Fit and the Best Fitting Model Using the Comparison to    ")
print("                   The Chi-Squared Distribution Using an alpha of 0.05                ")
print("======================================================================================")

#=================================================================================================|
#---Set up the matrix which will hold the goodness of fit results and then saved as a .csv file---|
#---And also set up the chi-squared critical values based on the degrees of freedom dof...--------|
#=================================================================================================|

goodness_fit_results <- matrix(ncol=5,nrow=2)
em = nrow(DR_Data); dofBP = em-2; dofExp = em-1;
gofstatBP = round(qchisq(0.95,dofBP),4); gofstatExp = round(qchisq(0.95,dofExp),4); 
Exp_goodfit_pvalue <- 1-pchisq(YE,dofExp)
BP_goodfit_pvalue <- 1-pchisq(YBP,dofBP)

#=====================================================================================|
#-------if loops in order to determine which conclusion to make, good fit or not------|
#----Then fill the matrix generated above with the necessary vales and conclusions----|
#=====================================================================================|
if (YBP<gofstatBP) {gofBPgood=("beta Poisson model shows a good fit to the data")} else{gofBPno=("beta Poisson model does not show a good fit to the data")}

if (YE<gofstatExp) {gofExpgood=("Exponential model shows a good fit to the data")} else{gofExpno=("Exponential model does not show a good fit to the data")}

goodness_fit_results[1,1] = ("Exponential"); goodness_fit_results[2,1] = ("Beta Poisson");
goodness_fit_results[1,2] = round(YE,4); goodness_fit_results[2,2] = round(YBP,4); goodness_fit_results[1,3] = gofstatExp
goodness_fit_results[2,3] = gofstatBP; goodness_fit_results[1,4] = round(Exp_goodfit_pvalue,4); goodness_fit_results[2,4] = round(BP_goodfit_pvalue,4)
if (YE<gofstatExp) {goodness_fit_results[1,5] = gofExpgood} else{goodness_fit_results[1,5] = gofExpno}
if (YBP<gofstatBP) {goodness_fit_results[2,5] = gofBPgood} else{goodness_fit_results[2,5] = gofBPno}
colnames(goodness_fit_results)<- c("MODEL","MINIMIZED DEVIANCE","CHI-SQUARED CRITICAL","CHI-SQRD P-value","CONCLUSION");
rownames(goodness_fit_results)<- c(" ", " ")

#==================================================================|
#----Make a .csv files which will be the goodness of fit matrix----|
#==================================================================|
write.csv(goodness_fit_results, file="Goodness_Fit_Results.csv")

#======================================================================|
#----Same as just previously but for determining best fitting model----|
#======================================================================|
best_fitting_model <- matrix(ncol=6,nrow=2)
bestmdlstat = round(qchisq(0.95,1),4)
deltaBPExp = abs(YE-YBP)
best_fit_pvalue <- 1-pchisq(deltaBPExp,1)

if (deltaBPExp > bestmdlstat) {bestfitbothBP=("beta Poisson model is the BEST fitting model")} else{bestfitbothExp=("Exponential is the BEST fitting model")}
best_fitting_model[1,1]=("Exponential"); best_fitting_model[2,1]=("Beta Poisson"); best_fitting_model[1,2]=round(YE,4)
best_fitting_model[2,2]=round(YE,4); best_fitting_model[2,2]=round(YBP,4); best_fitting_model[1,3]=round(deltaBPExp,4)
best_fitting_model[2,3]=round(deltaBPExp,4); best_fitting_model[1,4]=round(bestmdlstat,4); best_fitting_model[2,4]=round(bestmdlstat,4)
best_fitting_model[1,5]=round(best_fit_pvalue,4); best_fitting_model[2,5]=round(best_fit_pvalue,4)
if (deltaBPExp < bestmdlstat) {best_fitting_model[1,6]=bestfitbothExp} else{best_fitting_model[1,6]=bestfitbothBP}
if (deltaBPExp > bestmdlstat) {best_fitting_model[2,6]=bestfitbothBP} else{best_fitting_model[2,6]=bestfitbothExp}

colnames(best_fitting_model)<-c("MODEL","MINIMIZED DEVIANCE","DIFFERENCE BETWEEN DEVIANCES","CHI-SQUARED CRITICAL","CHI-SQRD P-value","CONCLUSION")
rownames(best_fitting_model)<-c(" "," ")
write.csv(best_fitting_model, file="Best_Fitting_Model.csv")
#if (YBP<gofstatBP) {print("yes")} else{print("no")}

#================================ Line 180 ================================|

#==========================================================================|
#---------Run bootstrap routine for number of dictated iterations----------|
#----------First is the exponential model then the beta Poisson------------|
#==========================================================================|
lk <- log(EXP_MLE[,2])
la <- log(BP_MLE[,2])
lN50 <- log(BP_MLE[,3])

iterations=10000
bootparms<-matrix(nrow=iterations,ncol=4)
for (iter in 1:iterations) 
	{
      bootdataframe=DR_Data
	total=bootdataframe$positive_response+bootdataframe$negative_response
	fobs=bootdataframe$positive_response/total
	bootpos<-rbinom(0+fobs,total,fobs)  # draw random sample
	bootdataframe$positive_response<-bootpos          # replace and form bootstrap sample

	results_boot<-mle(deviance.expntl,start=list(logk=lk), method = 'BFGS',
               fixed = list(obspos=bootdataframe$positive_response,
               obsneg=total-bootdataframe$positive_response,
               dose=bootdataframe$dose))
        results_boot
	L<-logLik(results_boot)
	L<-2*L[[1]]
	jb<-coef(results_boot)
      logk_est <- jb["logk"]
	k_est <- exp(jb["logk"])
	ExpID50=(log(-1/(0.5-1))/k_est)
	
	bootparms[iter,1] <- logk_est
	bootparms[iter,2] <- k_est 
	bootparms[iter,3] <- L
	bootparms[iter,4] <- ExpID50
 	}    

colnames(bootparms)<-c("ln(k)","k parameter","-2 ln(Likelihood)","LD50 or ID50")

n=10000
bootparms_bp<-matrix(nrow=n,ncol=6)
#bootbp <- matrix(nrow=n*6,ncol=1)
for (iter2 in 1:n) 
	{
      bootdataframe=DR_Data
	total=bootdataframe$positive_response+bootdataframe$negative_response
	fobs=bootdataframe$positive_response/total
	bootpos<-rbinom(0+fobs,total,fobs)  # draw random sample
	bootdataframe$positive_response<-bootpos          # replace and form bootstrap sample
	results_boot_bp<-mle(dev.bp,start=list(logalpha=la, logN50=lN50),
      method = 'BFGS', fixed = list(obspos=bootdataframe$positive_response,
      obsneg=total-bootdataframe$positive_response,dose=bootdataframe$dose))
      results_boot_bp
	LL<-logLik(results_boot_bp)
	LL<-2*L[[1]]
	jbp<-coef(results_boot_bp)
      logN50_est <- jbp["logN50"]
	logalpha_est <- jbp["logalpha"]
	N50_est <- exp(jbp["logN50"])
	alpha_est <- exp(jbp["logalpha"])
	BPID10=abs(((0.9^(-1/alpha_est)-1)*N50_est)/(2^(1/alpha_est)-1))
	BPID50=abs(((0.5^(-1/alpha_est)-1)*N50_est)/(2^(1/alpha_est)-1))	
	
	test <- bp.dr(alpha_est,N50_est,bootdataframe$dose)
	
	bootparms_bp[iter2,1] <- logalpha_est
	bootparms_bp[iter2,2] <- alpha_est
	bootparms_bp[iter2,3] <- logN50_est
	bootparms_bp[iter2,4] <- N50_est
	bootparms_bp[iter2,5] <- LL
	bootparms_bp[iter2,6] <- BPID10
#	bootbp[iter2,1] <- test
	}    

colnames(bootparms_bp)<-c("ln(alpha)","alpha","ln(N50)","N50","-2 ln(Likelihood)","LD10 or ID10")

#==================================Line 243==================================|

#============================================================================|
#-------Plot the results of the bootstrap routine of the exponential---------|
#============================================================================|

#==========histogram of exponential k parameters from the bootstrap==========|
	##pdf(paste("k_parameter_histogram.pdf",sep=''))
	png(paste("k_parameter_histogram.png",sep=''))
	par(cex=1, cex.axis=1.25,mai=c(1.0,1.2,0.58,0.15), cex.lab=1.5,font.lab=1)
	hist(bootparms[,1],breaks=20,plot=TRUE,xlab="ln(k)",main=" ",col="grey54")
	##dev.off()
	dev.off()
#==========Plot the exponential model with confidence intervals=============|

ndiv <- 700
klist <- exp(bootparms[,1]); write.csv(klist, file="klist.csv")
dmin <- min(dose)/100
ldmin <- log10(dmin)
dmax <- max(dose)*10
ldmax <- log10(dmax)
diff <- (ldmax - ldmin)/ndiv
bestfit<-matrix(nrow=ndiv+1,ncol=1)
plot01<-matrix(nrow=ndiv+1,ncol=1)
plot05<-matrix(nrow=ndiv+1,ncol=1)
plot95<-matrix(nrow=ndiv+1,ncol=1)
plot99<-matrix(nrow=ndiv+1,ncol=1)
plotdose <- matrix(nrow=ndiv+1,ncol=1)

for (iter in 0:ndiv+1) 
	{
   	pdose <- 10^(ldmin + (iter-1)*diff)
   	plotdose[iter] <- pdose
   	bestfit[iter] <- expntl.dr(exp(logk),pdose)
   	spread<- expntl.dr(klist, pdose)
   	CIs <- quantile(spread,probs=c(0.005,0.025,0.975,0.995)) 
   	plot01[iter] <- CIs[1]
   	plot05[iter] <- CIs[2]
   	plot95[iter] <- CIs[3]
   	plot99[iter] <- CIs[4]
	}

##pdf(paste("Exponential_model_curve.pdf",sep=''))
png(paste("Exponential_model_curve.png", sep=''))
par(cex=1,cex.lab=1.5,mai=c(1.2, 1.0,0.15,0.27),
    cex.axis=1.25,mgp=c(2.5,0.75,0), font.lab=1)
plot(dose,positive/(positive+negative),log="x",xlab="Dose (CFU)",
     ylab="Probability of Response",ylim=c(0,1),pch=17, xlim=c(dmin,dmax))

lines(plotdose,bestfit,lwd=2)
lines(plotdose,plot01,lty=2,lwd=2)
lines(plotdose,plot05,lty=3,lwd=2)
lines(plotdose,plot95,lty=3,lwd=2)
lines(plotdose,plot99,lty=2,lwd=2)
par(cex=1.3)
legend(dmin,1,legend=c("Exponential model","95% confidence","99% confidence"),
     lty=c(1,3,2), cex=0.7,bty="n",y.intersp=1.4)
##dev.off()
dev.off()
#=====================================Line 295=================================|

#==============================================================================|
#---------------Output the 95 and 99% CIs for the parameter (k)----------------|
#==============================================================================|
logk_CI <- quantile(bootparms[,1],probs=c(0.005,0.025,0.05,0.95,0.975,0.995))
k_CI <- data.frame(exp(logk_CI))
colnames(k_CI) <- ("k Confidence Interval")
write.csv(k_CI, file="k_Confidence_Intervals.csv")
EXPID50 <- quantile(bootparms[,4],probs=c(0.005,0.025,0.05,0.95,0.975,0.995))
EID50 <- data.frame(EXPID50)
colnames(EID50) <- ("ID50 Confidence Interval")
write.csv(EID50, file="Exponential_ID50_Confidence_Intervals.csv")

print("====================================================")
print("Confidence intervals for the exponential k parameter")
print(k_CI)
#============================================================================|

#============================================================================|
#-------Plot the results of the bootstrap routine of the beta Poisson--------|
#============================================================================|

#=========Scatter plot of beta Poisson parameters from the bootstrap=========|
#pdf(paste("alpha_N50_plot.pdf",sep=''))
#par(cex=1.2,mai=c(1.2,1.2,0.5,0.15),font.lab=1)
#plot(bootparms_bp[,1],bootparms_bp[,3],xlab=expression(ln(alpha)),ylab=expression(ln(N[50])))
#dataEllipse(bootparms_bp[,1],bootparms_bp[,3],levels=c(0.9,0.95,0.99), center.pch=0,col="black",lwd=0.5, lty=c(1,2,3),plot.points=F, add=T)
#dev.off()

##pdf(paste("alpha_N50_plot.pdf",sep=''))
png(paste("alpha_N50_plot.png",sep=''))
##par(cex=1.2,mai=c(1.2,1.2,0.5,0.15),font.lab=1)
par(cex=1,cex.lab=1.5,mai=c(0.8,0.8,0.4,0.2),cex.axis=1.25,mgp=c(2.5,0.75,0),font.lab=1)
plot(bootparms_bp[,1],bootparms_bp[,3],xlab=expression(ln(alpha)),ylab=expression(ln(N[50])))
dataEllipse(bootparms_bp[,1],bootparms_bp[,3],levels=0.9, center.pch=0,col="red",lwd=2, lty=1,plot.points=F, add=T)
dataEllipse(bootparms_bp[,1],bootparms_bp[,3],levels=0.95, center.pch=0,col="green",lwd=2, lty=2,plot.points=F, add=T)
dataEllipse(bootparms_bp[,1],bootparms_bp[,3],levels=0.99, center.pch=0,col="blue",lwd=2, lty=4,plot.points=F, add=T)
points(x=logalpha, y=logN50, pch=21, cex=1.5, bg="gray50")
points(x=logalpha, y=logN50, pch=4, cex=3, font=2)
legend("top",legend=c("0.90 Confidence", "0.95 Confidence", "0.99 Confidence"),cex=0.90,horiz="TRUE",col=c("red","green","blue"),lwd=c(1.5,1.5,1.5), lty=c(1,2,4),
pt.cex=2, x.intersp=0.1, xjust=1, box.col=NULL, bg="white", bty='n', xpd=TRUE, inset=-0.07)

##legend("top",legend=c("Parameter MLEs", "0.90 Confidence", "0.95 Confidence", "0.99 Confidence"),cex=0.90,horiz="TRUE", 
##lty=c(0,1,2,3), pch=c(13, NA, NA, NA), pt.cex=2, x.intersp=0.1, xjust=1, box.col=NULL, bg="white", bty='n', xpd=TRUE, inset=-0.07) 
##legend("top",legend=c("0.90 Confidence", "0.95 Confidence", "0.99 Confidence"),cex=.75,horiz="TRUE", lty=c(1,2,4),bg="white")
##dev.off()
dev.off()

#===========Plot the beta Poisson model with confidence intervals============|

ndiv <- 700
alphalist <- exp(bootparms_bp[,1]); write.csv(alphalist, file = "alphalist.csv")
N50list <- exp(bootparms_bp[,3]); write.csv(N50list, file="N50list.csv")
dmin <- min(dose)/100
ldmin <- log10(dmin)
dmax <- max(dose)*10
ldmax <- log10(dmax)
diff <- (ldmax - ldmin)/ndiv
bestfit_bp<-matrix(nrow=ndiv+1,ncol=1)
plot01_2<-matrix(nrow=ndiv+1,ncol=1)
plot05_2<-matrix(nrow=ndiv+1,ncol=1)
plot95_2<-matrix(nrow=ndiv+1,ncol=1)
plot99_2<-matrix(nrow=ndiv+1,ncol=1)
plotdose <- matrix(nrow=ndiv+1,ncol=1)

for (iter in 0:ndiv+1) {
   pdose <- 10^(ldmin + (iter-1)*diff)
   plotdose[iter] <- pdose
   bestfit_bp[iter] <- bp.dr(exp(logalpha),exp(logN50),pdose)
   spread <- bp.dr(alphalist,N50list, pdose)
   CI <- quantile(spread,probs=c(0.005,0.025,0.975,0.995)) 
   plot01_2[iter] <- CI[1]
   plot05_2[iter] <- CI[2]
   plot95_2[iter] <- CI[3]
   plot99_2[iter] <- CI[4]}

##pdf(paste("Beta_Poisson_model_curve.pdf",sep=''))
png(paste("Beta_Poisson_model_curve.png", sep=''))
par(cex=1.0,cex.lab=1.4,mai=c(1.2, 1.2,0.15,0.3),
    cex.axis=1.25,mgp=c(2.5,0.75,0), font.lab=1)
plot(dose,positive/(positive+negative),log="x",xlab=expression(Dose (CFU)),
     ylab="Probability of Response",ylim=c(0,1),pch=17, xlim=c(dmin,dmax))

lines(plotdose,bestfit_bp,lwd=2)
lines(plotdose,plot01_2,lty=2,lwd=2)
lines(plotdose,plot05_2,lty=3,lwd=2)
lines(plotdose,plot95_2,lty=3,lwd=2)
lines(plotdose,plot99_2,lty=2,lwd=2)
par(cex=1.3)
legend(dmin,1,legend=c("beta Poisson model","95% confidence","99% confidence"),
     lty=c(1,3,2),cex=0.8,bty="n",y.intersp=1.4)
##dev.off()
dev.off()
#====================================Line 361====================================|

#================================================================================|
#-----------------Output the 95 and 99% CIs for alpha and N50--------------------|
#================================================================================|

alpha_CI <- quantile(bootparms_bp[,1],probs=c(0.005,0.025,0.05,0.95,0.975,0.995))
ealpha_CI <- data.frame(exp(alpha_CI))
colnames(ealpha_CI)<-("alpha confidence Interval")
write.csv(ealpha_CI, file="alpha_Confidence_Intervals.csv")
print("=========================================================")
print("Confidence intervals for the beta Poisson alpha parameter")
print(ealpha_CI)

N50_CI <- quantile(bootparms_bp[,3],probs=c(0.005,0.025,0.05,0.95,0.975,0.995))
eN50_CI <- data.frame(exp(N50_CI))
colnames(eN50_CI)<-("N50_Confidence_Intervals")
write.csv(eN50_CI, file="N50_Confidence_Interval.csv")
print("=========================================================")
print("Confidence intervals for the beta Poisson N50 parameter")
print(eN50_CI)

BPID50_CI <- quantile(bootparms_bp[,6], probs=c(0.005,0.025,0.05,0.95,0.975,0.995))
dBPID50_CI <- data.frame(BPID50_CI)
colnames(dBPID50_CI)<-("LD50 or ID50")
write.csv(dBPID50_CI, file="LD_or_ID50_betaPoisson.csv")
print(dBPID50_CI)

}else{print("There is NOT a trend between dose and observed probability of response")}#end of if statement to perform the trend test
#=================================END OF CODE=====================================|
#dev.off()