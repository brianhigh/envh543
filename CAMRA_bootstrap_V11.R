#==========================================================================|
#-----------------------------------HEADER---------------------------------|
#==========================================================================|
#CAMRA_bootstrap_v11.r is an R statistical programming package		   |
#source code for initializing the bootstrap technique in R for		   |
#generation of confidence intervals for both the exponential and           |
#beta Poisson dose response models for risk estimation of                  |
#exposure to pathogens                                                     |
#Written for CAMRA by: Mark H. Weir, Timothy A. Bartrand and               |
#Charles N. Haas, Department of the Civil Architectural and Environmental  |
#Engineering, Drexel University                                            |
#==========================================================================|

# == V11 Update ===
# This version uses the gsl package to allow for the hypergeometric        | 
# function to be used. If you are on a mac you will need to use HomeBrew   |
# to install the gsl package before running this code. If you are on a PC  |
# or Linux this is not needed for those standard distros tests - tested on |
# Fedora 32+; Ubuntu 18+; and  SUSE					   |
#==========================================================================|						

#==========================================================================|
#--Draw in dose response data and assign the needed values from the data---|
#==========================================================================|
#DR_Data <- read.table("BA_RM_NEW_format.txt",header=TRUE)				  
  
  dose <- DR_Data$dose									  
	positive <- DR_Data$positive_response									  
	negative <- DR_Data$negative_response									   
	oprob <- positive/(positive+negative)			
	BP_Iter <- factorial(length(DR_Data$dose))
#==========================================================================|

#==========================================================================|
#-----------Load required libraries to perform MLE and bootstrap-----------|
#==========================================================================|
	require(stats4)									         
	require(boot)	
	require(car)	
	require(gsl)
#==========================================================================|

dose = DR_Data$dose
positive = DR_Data$positive_response
negative = DR_Data$negative_response
ni = positive+negative

xi = log(dose)
xbar = sum(ni*xi)/sum(ni)
pbar = sum(positive)/sum(ni)


#==========================================================================|
#----------------------Define the dose response models.--------------------|
#----expntl.dr ---> exponential model and bp.dr ---> beta Poisson model----|
#==========================================================================|
	expntl.dr <- function(k,dose) 1 - exp(-k*dose)    
	bp.dr <- function(alpha,N50,dose) 1-(1+(dose/N50)*(2^(1/alpha)-1))^(-alpha)
	ex.bp.dr <- function(a,b,dose) 1 - hyperg_1F1(a,a+b,-dose)
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

dev.ex.bp <- function (obspos, obsneg, loga, logb, dose)
  {
    eps = 1e-15;
    a = exp(loga)
    b = exp(logb)
    obsf = obspos/(obspos+obsneg);
    pred = ex.bp.dr(a,b,dose);
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
	write.csv(EXP_MLE, file="Exponential_MLE_Output.csv")
	

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
	write.csv(BP_MLE, file = "betaPoisson_MLE_Output.csv", row.names = FALSE)
	
results_exbp <- mle(dev.ex.bp,start=list(loga=-0.2, logb=10),
                    method='BFGS',fixed=list(obspos=positive,obsneg=negative,dose=dose))
  EX_BP_MLE <- matrix(ncol=3,nrow=1)
  Y_EX_BP <- logLik(results_exbp)
  EX_BP_LL <- -2*Y_EX_BP[[1]]
  jc <- coef(results_exbp)
  EX_BP_alpha <- exp(jc["loga"])
  EX_BP_beta <- exp(jc["logb"])

  EX_BP_MLE[,1] <- EX_BP_LL
  EX_BP_MLE[,2] <- EX_BP_alpha
  EX_BP_MLE[,3] <- EX_BP_beta

  colnames(EX_BP_MLE)<- c("Minimized Deviance","alpha","beta"); 
  rownames(EX_BP_MLE)<-(" ")
  write.csv(EX_BP_MLE, file = "EXACT_betaPoisson_MLE_Output.csv", row.names = FALSE)

	
#=================================================================================================|
#---Set up the matrix which will hold the goodness of fit results and then saved as a .csv file---|
#---And also set up the chi-squared critical values based on the degrees of freedom dof...--------|
#=================================================================================================|

goodness_fit_results <- matrix(ncol=5,nrow=2) 
em = nrow(DR_Data); dofBP = em-2; dofExp = em-1;
gofstatBP = round(qchisq(0.95,dofBP),4); gofstatExp = round(qchisq(0.95,dofExp),4); 
Exp_goodfit_pvalue <- 1-pchisq(YE,dofExp)
BP_goodfit_pvalue <- 1-pchisq(EX_BP_LL,dofBP)

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
write.csv(goodness_fit_results, file = "Goodness_Fit_Results.csv", row.names = FALSE)

#==================================================================|
#----Make a .csv files which will be the goodness of fit matrix----|
#==================================================================|

#======================================================================|
#----Same as just previously but for determining best fitting model----|
#======================================================================|
best_fitting_model <- matrix(ncol=6,nrow=2)
bestmdlstat = round(qchisq(0.95,1),4)
deltaBPExp = abs(YE-EX_BP_LL)
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
write.csv(best_fitting_model, file = "Best_Fitting_Model.csv", row.names = FALSE)


#==================================================================|
#---- This is where the bootstrap simulation would be          ----|
#==================================================================|

n=BP_Iter 
bootparms_bp<-matrix(nrow=n,ncol=6)
for (iter2 in 1:n) 
{
  bootdataframe=DR_Data
  total=bootdataframe$positive_response+bootdataframe$negative_response
  fobs=bootdataframe$positive_response/total
  bootpos<-rbinom(0+fobs,total,fobs)  # draw random sample
  bootdataframe$positive_response<-bootpos          # replace and form bootstrap sample
  results_boot_bp<-mle(dev.bp,start=list(logalpha=0.1, logN50=0.5),method='L-BFGS-B',lower=rep(-2,2),upper=rep(3.5,2),
                       control=list(trace=0),
                       fixed = list(obspos=bootdataframe$positive_response,
                                    obsneg=total-bootdataframe$positive_response,dose=bootdataframe$dose))
  results_boot_bp
  LL<-logLik(results_boot_bp)
  LL<-2*LL[[1]]
  jbp<-coef(results_boot_bp)
  logN50_est <- jbp["logN50"]
  logalpha_est <- jbp["logalpha"]
  N50_est <- exp(jbp["logN50"])
  alpha_est <- exp(jbp["logalpha"])
  BPID10=abs(((0.9^(-1/alpha_est)-1)*N50_est)/(2^(1/alpha_est)-1))
  
  bootparms_bp[iter2,1] <- logalpha_est
  bootparms_bp[iter2,2] <- alpha_est
  bootparms_bp[iter2,3] <- logN50_est
  bootparms_bp[iter2,4] <- N50_est
  bootparms_bp[iter2,5] <- LL
  bootparms_bp[iter2,6] <- BPID10
}    

colnames(bootparms_bp)<-c("ln(alpha)","alpha","ln(N50)","N50","-2 ln(Likelihood)","LD50 or ID50")
write.csv(bootparms_bp, file="Approx_BetaPoisson_Bootstrap.csv")

n=BP_Iter 
bootparms_EX_BP<-matrix(nrow=n,ncol=5)
for (iter2 in 1:n) 
{
  bootdataframe=DR_Data
  total=bootdataframe$positive_response+bootdataframe$negative_response
  fobs=bootdataframe$positive_response/total
  bootpos<-rbinom(0+fobs,total,fobs)  # draw random sample
  bootdataframe$positive_response<-bootpos          # replace and form bootstrap sample
  results_boot_EXbp<-mle(dev.ex.bp,start=list(loga=0.1, logb=0.5),method='BFGS',
                       control=list(trace=0), fixed = list(obspos=bootdataframe$positive_response,
                                    obsneg=total-bootdataframe$positive_response,dose=bootdataframe$dose))
  results_boot_EXbp
  LL<-logLik(results_boot_EXbp)
  LL<-2*LL[[1]]
  jbp<-coef(results_boot_EXbp)
  logb_est <- jbp["logb"]
  loga_est <- jbp["loga"]
  b_est <- exp(jbp["logb"])
  a_est <- exp(jbp["loga"])

  bootparms_EX_BP[iter2,1] <- loga_est
  bootparms_EX_BP[iter2,2] <- a_est
  bootparms_EX_BP[iter2,3] <- logb_est
  bootparms_EX_BP[iter2,4] <- b_est
  bootparms_EX_BP[iter2,5] <- LL
}    

colnames(bootparms_EX_BP)<-c("ln(alpha)","alpha","ln(Beta)","Beta","-2 ln(Likelihood)")
write.csv(bootparms_EX_BP, file="Exact_BetaPoisson_Bootstrap.csv")


