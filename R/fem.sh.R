fem.computeSH <- function(out,plot=TRUE){
  # Compute SH criterion
  dd = data.frame(nbp=out$allCriteria$comp,ll=out$allCriteria$loglik)
  fit = rlm(ll  ~ nbp,data=dd,method='MM')
  if(fit$coefficients[2]<0) fit$coefficients[2]=0
  dd$llpen = dd$ll- 2* fit$coefficients[2]*dd$nbp
  list(opt=which.max(dd$llpen),crit=dd$llpen)
  
  # Plot
  old.par <- par(no.readonly = TRUE)
  par(mfrow=c(1,2))
  plot(dd$nbp,dd$ll,type='p',xlab='Model complexity',ylab='Log-likelihood')
  abline(fit,col='red')
  title('Slope computation')
  plot(dd$nbp,dd$llpen,type='p',xlab='Model complexity',ylab='SH criterion')
  points(dd$nbp[which.max(dd$llpen)],max(dd$llpen),col='red',pch=19)
  text(dd$nbp[which.max(dd$llpen)],max(dd$llpen),
  labels=paste(as.character(out$allCriteria$model[which.max(dd$llpen)]),' K=',
                out$allCriteria$K[which.max(dd$llpen)],sep=''),pos = 4,cex=0.8)
  title('Slope heuristic criterion')
  par(old.par)
  
  # Select best model according to SH
  id_max = which.max(dd$llpen)
  if (is.matrix(out$allResults)) {res = out$allResults[,id_max]} else res = out$allResults[[id_max]]
  res$call = out$call; class(res) = 'fem'
  return(res)
}