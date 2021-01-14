library(survival)
library(matlabr)
pbc.mydata = pbc[is.na(pbc$protime) == F,]
pbc.mydata = pbc.mydata[is.na(pbc.mydata$albumin) == F,]
pbc.mydata = pbc.mydata[is.na(pbc.mydata$bili) == F,]
pbc.mydata = pbc.mydata[is.na(pbc.mydata$age) == F,]
pbc.mydata = pbc.mydata[is.na(pbc.mydata$edema) == F,]

# Changing the time to year scale
pbc.mydata$time_yr <- pbc.mydata$time/365.25

##########################################################################################
# You should construct data matrix(i.e. m_data): first column should be log(survival time) 
#                                               last column should be censoring indicator
##########################################################################################
m_data=cbind(log(pbc.mydata$time_yr), pbc.mydata$age,pbc.mydata$edema,log(pbc.mydata$bili),
             log(pbc.mydata$protime),log(pbc.mydata$albumin),pbc.mydata$status==2)

# The below is to call matlab function, save the result, and read the result.
# You do not have to modify the below. The result will be stored in MAFT
write.table(m_data, file="./tmp_data.csv", sep=",", row.names=FALSE, col.names=FALSE)
run_matlab_script("matlab_fit.m")
coefficient=as.matrix(read.csv("./tmp_result.csv", header=FALSE)) # estimated beta
mixing=as.matrix(read.csv("./tmp_result_mixing.csv", header=FALSE)) 
cons=as.numeric(read.csv("./tmp_result_cons.csv", header=FALSE))
rownames(mixing)=c('support','weight')
colnames(mixing)=NULL
result=list(coefficients=coefficient,mixing=mixing,cons=cons)
##########################################################################
#result$coefficients: beta estimator
#result$mixing: support and weight for the estimated mixing distribution
#result$cons: the lower bound of the support for the mixing distribution
