# library(boostnew)
# library(glmnet)
# library(pROC)
# library(glmnet)
# library(parallel)
# library("MASS")
# library("Matrix")
# require(compiler)
# library(IGESS)
# library(Hmisc)
# library(iGESSSimu)
# enableJIT(4)
# # library(boostnew)
# # library(Matrix)
# # # crossvalind <- function(N, kfold) {
# # #   len.seg <- ceiling(N/kfold)
# # #   incomplete <- kfold*len.seg - N
# # #   complete <- kfold - incomplete
# # #   ind <- matrix(c(sample(1:N), rep(NA, incomplete)), nrow = len.seg, byrow =
# # #                   TRUE)
# # #   cvi <- lapply(as.data.frame(ind), function(x) c(na.omit(x))) # a list
# # #   return(cvi)
# # # }
# # #protaste_cancer
# # X <- read.csv("/Users/david_dai/Downloads/real data/protaste_cancer/X.csv",header = F)
# # y <- read.csv("/Users/david_dai/Downloads/real data/protaste_cancer/y.csv",header = F)
# # X <- as.matrix(X)
# # #y <- as.integer(y)
# # y <- as.matrix(y)
# #
# # fit0 <- boostnew::LassoPathCV4Interact(X, y, 0)
# # fit100 <- boostnew::LassoPathCV4Interact(X, y, 100)
# # min(fit0$cvm)
# # min(fit100$cvm)
# # yhat <- X %*% fit0$beta[,fit0$index_min] + fit0$covs[fit0$index_min]
# #
# # plotPath(fit100)
# # plotCVFit(fit0)
#
# #Leukemia
# data <- read.csv("/Users/david_dai/Downloads/real data/Leukemia/leukemia_small.csv",header = F)
# data <- t(data)
# y <- matrix(0,nrow(data),1)
# y[which(data[,1] == "ALL")] = 1
# y[which(data[,1] == "AML")] = 0
# X <- as.matrix(data[,-1])
# X <- X[-1,]
# y <- y[-1]
# N = nrow(X)
# P = ncol(X)
# X<-matrix(as.numeric(X),N,P)
# L<-BOLTSSIRR::BOLT_SSI_RR(X, y)
#


data <- read.csv("/Users/david_dai/Downloads/real data/supermarket/market.csv",header = F)
X <- data[,-1]
X <- as.matrix(X)
y <- data[,1]
y <- as.matrix(y)
library(BOLTSSIRR)
L <- CV_BOLT_SSI_RR(X,y)

setwd("BioBank")
text <- readLines("/Users/david_dai/Desktop/1419.txt")
L <-strsplit(text, "</a></td>")
path <- "https://s3.amazonaws.com/gene2pheno/"
idx <- 0
for(i in 9:22){
  str = L[[i]]
  subL <-strsplit(str, "</a></td>")
  print(length(subL))
  for(j in 1:1){

    subsubL <-strsplit(subL[[j]], "\">")
    X <- subsubL[[1]]
    substr <- X[[4]]
    idx <- idx + 1
    print(paste(substr,idx))
    # if(i== 10 && j == 17){
    #   print(substr);
    # }
   # if(i > 11 || j >=43)
    # if(file.exists(substr) == FALSE)
    # {
    #   cmd <- (paste0("wget ",path,substr))
    #   system(cmd)
    #   print(paste(idx, substr,"created!"));
    #
    # }else{
    #   print(paste(idx, substr,"existed!"));
    #   idx <- idx + 1
    # }



  }

}



load("interaction_NFBC.Rdata")
L <- BOLTSSIRR::BOLT_SSI(x[,1:1000], y, extra_pairs = 5000, thread_num = 4)

load("interaction_NFBC.Rdata")
L <- BOLTSSIRR::BOLT_SSI(x[,1:100], y, extra_pairs = 5000, thread_num = 4)
#L[[8]] to L[[22]]
L <- BOLTSSIRR::BOLT_SSI(x, y, extra_pairs = 10000, thread_num = 22)
