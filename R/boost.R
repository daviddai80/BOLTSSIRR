plotPath <- function(fit){
  which=nonzeroCoef(fit$beta)
  nwhich=length(which)
  beta=as.matrix(fit$beta[which,,drop=FALSE])
  index=log(fit$lambdas)
  iname="Log Lambda"
  approx.f=0
  xlab = iname
  ylab="Coefficients"
  l <- log(fit$lambdas)
  matplot(index,t(beta),lty=1,type="l",xlab = xlab, ylab=ylab)#, xlim=rev(range(l))
 # max(fit$lambdas)
  z = colSums(beta != 0)
  atdf=pretty(index)
  prettydf=approx(x=index,y=z,xout=atdf,rule=2,method="constant",f=approx.f)$y
  axis(3, at = atdf,labels = prettydf,tcl=NA)
}

plotCVFit <- function(fit){
   cvobj = fit;
   sign.lambda = 1
   which=nonzeroCoef(cvobj$beta)
   beta=as.matrix(cvobj$beta[which,,drop=FALSE])
   nonzero = colSums(beta != 0)
   xlab="log(Lambda)"
   if(sign.lambda<0)xlab=paste("-",xlab,sep="")
   plot.args=list(x=sign.lambda*log(cvobj$lambdas),y=cvobj$cvm,ylim=range(cvobj$cvup,cvobj$cvlo),xlab=xlab,ylab=cvobj$name,type="n")
   # new.args=list(...)
   # if(length(new.args))plot.args[names(new.args)]=new.args
   do.call("plot",plot.args)
   errbar(sign.lambda*log(cvobj$lambdas), cvobj$cvm,cvobj$cvup,cvobj$cvlo)
   points(sign.lambda*log(cvobj$lambdas),cvobj$cvm,pch=20,col="red")
   axis(side=3,at=sign.lambda*log(cvobj$lambdas),labels=paste(nonzero),tick=FALSE,line=0)
   abline(v=sign.lambda*log(cvobj$lambda_min),lty=3)
   abline(v=sign.lambda*log(cvobj$lambda_1se),lty=3)
   invisible()

}
