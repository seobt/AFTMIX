library(survival)
library(aftgee)
library(rms) # for bj
library(matlabr)
data(pbc)
head(pbc)
pbc.mydata = pbc[is.na(pbc$protime) == F,]
pbc.mydata = pbc.mydata[is.na(pbc.mydata$albumin) == F,]
pbc.mydata = pbc.mydata[is.na(pbc.mydata$bili) == F,]
pbc.mydata = pbc.mydata[is.na(pbc.mydata$age) == F,]
pbc.mydata = pbc.mydata[is.na(pbc.mydata$edema) == F,]

# Changing the time to year scale
pbc.mydata$time_yr <- pbc.mydata$time/365.25

# Fitting a smoothed Gehan estimator: using time in years
SAFT_R = aftsrr(Surv(time_yr,status==2) ~ age + edema + 
                            log(pbc.mydata$bili) + log(pbc.mydata$protime) + 
                            log(pbc.mydata$albumin),data=pbc.mydata,B=0,se = "ISMB") #Gehan
######## Calculation of the intercept
data.x=as.matrix(cbind(pbc.mydata$age,pbc.mydata$edema,log(pbc.mydata$bili),
                       log(pbc.mydata$protime),log(pbc.mydata$albumin)))
mydata1 = log(pbc.mydata$time_yr) - data.x%*%as.matrix(SAFT_R$beta)
ind=ifelse(pbc.mydata$status==2,1,0)
mydata.surv <- survfit(Surv(mydata1, ind)~ 1, conf.type="none")
est.int = mydata.surv$time%*%diff(c(0,(1-mydata.surv$surv))) # Estimated intercept for Gehan

# Least square estimator using GEE
SAFT_G=aftgee(Surv(time_yr,status==2)~age + edema + 
              log(pbc.mydata$bili) + log(pbc.mydata$protime) + 
              log(pbc.mydata$albumin),data=pbc.mydata,B=0)

## Parametric MLE under normality
PAFT=survreg(Surv(time_yr,status==2)~age + edema + 
          log(pbc.mydata$bili) + log(pbc.mydata$protime) + 
          log(pbc.mydata$albumin),data=pbc.mydata,dist="lognormal")

##################################################################################
### Run Matlab code for one sample
#################################################################################
# You should construct data matrix(i.e. m_data): first column should be log(survival time) 
#                                               last column should be censoring indicator

m_data=cbind(log(pbc.mydata$time_yr), pbc.mydata$age,pbc.mydata$edema,log(pbc.mydata$bili),
             log(pbc.mydata$protime),log(pbc.mydata$albumin),pbc.mydata$status==2)

# You do not have to modify the below. The result will be stored in MAFT
# The first row of MAFT is estimated regression coefficients
write.table(m_data, file="./tmp_data.csv", sep=",", row.names=FALSE, col.names=FALSE)
run_matlab_script("matlab_fit.m")
MAFT <- as.matrix(read.csv("./tmp_result.csv", header=FALSE)) # estimated beta
mixing=as.matrix(read.csv("./tmp_result_mixing.csv", header=FALSE))
cons=as.numeric(read.csv("./tmp_result_cons.csv", header=FALSE))

# Compute P(survival time>t) given covariate
compute_prob=function(t,est_beta,covariate,est_mixing,est_cons)
{
  #est_beta: estimated beta
  #est_mixing: estimated mixing distribution 
  K=dim(est_mixing)[2]
  xx=cbind(1,covariate)
  prob=0
  for (j in 1:K){
    prob=prob+est_mixing[2,j]*(1-pnorm(t,xx%*%t(est_beta),sqrt(est_mixing[1,j]+est_cons)))
  }
  return(prob)
}

# The following example returns the estimated probability that the log(survival) time 
# is greater than 1 when the covariate is the same as the first observation in PBC data
compute_prob(1,MAFT,m_data[1,2:6],mixing,cons)



##################################################################################
### Run Matlab code for multiple samples
#################################################################################
mydata=read.csv(file("./normal200_5.csv"),header=F)
mydata=mydata[,1:20]
# mydata is the matrix that stacks each sample to right 
# For example, see "normal200_5.csv" which was used for simulation 1 in our manuscript
# The result is stored in "betas" in which each row represents beta estiamtes of each sample
write.table(mydata, file="./tmp_mydata.csv", sep=",", row.names=FALSE, col.names=FALSE)
run_matlab_script("matlab_sim.m")
betas <- as.matrix(read.csv("./tmp_beta.csv", header=FALSE)) #store estimated beta

