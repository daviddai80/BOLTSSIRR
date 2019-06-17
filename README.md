# BOLTSSIRR
BOLTSSI is a statistical approach for detecting interaction effects among predict variables to response variables, which is often an crucial step in regression modeling of real data for various applications. Through this publicly available package, we provide a unified environment to carry out interaction pursuit using a simple sure screening procedure (SSI) to fully detect significant pure interactions between predict variables and the response variable in the high or ultra-high dimensional generalized linear regression models. Furthermore, we suggest to discretize continuous predict variables, and utilize the Boolean operation for the marginal likelihood estimates. The so-called ‘BOLTSSI’ procedure is proposed to accelerate the sure screening speed of the procedure.



Usage
=======

The file Package_BOLTSSI.pdf in inst/doc shows several examplex of how to use BOLTSSIRR package. 

1. Sure Independence Screening SSI and BOLTSSI for interaction screening.
```R
library(BOLTSSIRR)
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

