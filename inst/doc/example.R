
rm(list=ls())
library(BOLTSSIRR)

set.seed(0)
p=300;n=100
####################generate the covariance matrix
rho=0.5
H<-abs(outer(1:p,1:p,"-"))
covxx=rho^H
cholcov = chol(covxx)

x0 = matrix(rnorm(n*p), n, p)
x = x0%*%cholcov
#gaussian response
set.seed(0)
y=2*x[,1]+2*x[,8]+3*x[,1]*x[,8]+rnorm(n)



model1=BOLT_SSI(x,y)
head(model1)


#binary response
set.seed(40)
feta = 2*x[,1]+2*x[,8]+3*x[,1]*x[,8]; 
fprob = exp(feta)/(1+exp(feta))
y = rbinom(n, 1, fprob)
model2=BOLT_SSI(x,y)
head(model2)



###########################
set.seed(0)
p=300;n=100
####################generate the covariance matrix
rho=0.5
H<-abs(outer(1:p,1:p,"-"))
covxx=rho^H
cholcov = chol(covxx)

x0 = matrix(rnorm(n*p), n, p)
x = x0%*%cholcov
#gaussian response
set.seed(0)
y=2*x[,1]+2*x[,8]+3*x[,1]*x[,8]+rnorm(n)
model3=CV_BOLT_SSI_RR(x,y,extra_pairs=p,nfold=5)


Lambdas=model3$lambdas
Beta=model3$beta
Lambda.min=model3$lambda_min
index.min=which(Lambdas==Lambda.min)
final_Beta=Beta[,index.min]
Pairs=t(matrix(model3$pairs,nrow = 2))

m=length(final_Beta)
main_index=which(final_Beta[1:p]!=0)
inter_index=which(final_Beta[(p+1):m]!=0)
final_inter_index=Pairs[inter_index,]

newx=matrix(rnorm(20*p),20)
predictY=BOLT_Predict(newx,model3)

