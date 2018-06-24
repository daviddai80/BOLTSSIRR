# BOLT_Predict <- function(X, fit){
#     v <- as.vector(fit$pairs)
#     npair <- nrow(fit$pairs) / 2
#     Xc <- matrix(0, nrow(X), npair)
#     for(i in 1:npair){
#        c1 <- v[((i-1)*2 + 1)]
#        c2 <- v[i*2]
#        Xc[,i] <- X[,c1] * X[,c2];
#     }
#     X <- cbind(X,Xc)
#     yhat = fit$cov[fit$index_min] + X * fit$beta[ , fit$index_min];
# }


BOLT_Predict <- function(X, fit){
  v <- as.vector(fit$pairs)
  V <- matrix(v,2)
  Xc <-X[,V[1,]]*X[,V[2,]]
  X <- cbind(X,Xc)
  yhat = fit$covs[fit$index_min] + X %*% fit$beta[ , fit$index_min];
}

