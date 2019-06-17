# BOLTSSIRR
BOLTSSI is a statistical approach for detecting interaction effects among predict variables to response variables, which is often an crucial step in regression modeling of real data for various applications. Through this publicly available package, we provide a unified environment to carry out interaction pursuit using a simple sure screening procedure (SSI) to fully detect significant pure interactions between predict variables and the response variable in the high or ultra-high dimensional generalized linear regression models. Furthermore, we suggest to discretize continuous predict variables, and utilize the Boolean operation for the marginal likelihood estimates. The so-called ‘BOLTSSI’ procedure is proposed to accelerate the sure screening speed of the procedure.



Usage
=======

The file Package_BOLTSSI.pdf in inst/doc shows several examplex of how to use BOLTSSIRR package. 

1. Sure Independence Screening SSI and BOLTSSI for interaction screening.
```R
library(BOLTSSIRR)
library(BOLTSSIRR) 
set.seed(0) 
p=300;
n=100;
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
```

Development 
=======
This R package is developed by Zhoumin, Mingwei Dai, Can Yang, Heng Peng, and maintained by Can Yang <eeyangc@gmail.com>.

Installation
=======
To install the development version of BOLTSSIRR, it's easiest to use the 'devtools' package. Note that BOLTSSIRR depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

install.packages("devtools")  
library(devtools)  
install_github("daviddaigithub/BOLTSSIRR")  

References
=======
Min Zhou, Mingwei Dai, Yuan Yao, Jin Liu, Can Yang, Heng Peng. BOLT-SSI: A Statistical Approach to Screening Interaction Effects for Ultra-High Dimensional Data. https://arxiv.org/abs/1902.03525

