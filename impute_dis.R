impute_dis=function(X,score){
  k=floor(sqrt(dim(score)[1]))
  sim=as.matrix(dist(as.matrix(score), method = "euclidean"))
  sim=max(sim)-sim;
  KK = colSums(sim);
  twom = sum(KK);
  sim_norm=sim - as.matrix(KK)%*%t(as.matrix(KK))/twom
  f=apply(sim_norm, 2, order);
  i=rep(1:dim(sim)[2],each = k)
  p=as.vector(f[1:k,])
  KNN=Matrix::sparseMatrix(i,p,x=rep(1/k,length(i)),dims=c(dim(sim)[2],dim(sim)[2]))
  impute=KNN%*%Matrix::t(X);
  return(impute)
}