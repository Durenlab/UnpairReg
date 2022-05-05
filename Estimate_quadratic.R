Estimate_quadratic=function(X,Y,beta,beta_nonzero_idx,beta_dis,d0,lambda){
  n=dim(Y)[1]
  Xbeta=matrix(0,nrow=dim(X)[1],ncol=n)
  for (i in 1:n){
    Xbeta[,i]=X[,beta_nonzero_idx[,i]==1]%*%beta[beta_nonzero_idx[,i]==1,i];
  }
  for (i in 1:n){
    TFAS=exp(beta_dis[beta_nonzero_idx[,i]==1,i]/d0)
    TFAS[1]=0
    W=as.matrix(TFAS)%*%t(as.matrix(TFAS))
    Temp=t(Xbeta[beta_nonzero_idx[,i]==1,setdiff(1:n,i)])
    beta[beta_nonzero_idx[,i]==1,i]=pinv(t(Temp)%*%Temp+lambda*W)%*%t(Temp)%*%Y[setdiff(1:n,i),i]#pinv inverse
    Xbeta[,i]=X[,beta_nonzero_idx[,i]==1]%*%beta[beta_nonzero_idx[,i]==1,i];
  }
  return(beta)
}
