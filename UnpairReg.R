UnpairReg=function(result,d0,lambda){
library(tidyr)
library(Matrix)
library(pracma)
E=result[[1]]
O=result[[2]]
Symbol_location=result[[3]]
Peak_location=result[[4]]
gene_count=Matrix::rowSums(E>0)
region_count=Matrix::rowSums(O>0)
N1=floor(dim(E)[2]*0.1);
N2=floor(dim(O)[2]*0.1);
TG_filter=E[gene_count>N1,]
RE_filter=O[region_count>N2,]
Symbol_location=Symbol_location[gene_count>N1,]
Peak_location=Peak_location[region_count>N2,]
TG_filter=TG_filter[,Matrix::colSums(TG_filter)>0]
RE_filter=RE_filter[,Matrix::colSums(RE_filter)>0]
A=crossing(1:dim(Symbol_location)[1],1:dim(Peak_location)[1])
colnames(A)=c('symbol','peak')
x=Symbol_location[A$symbol,1]
y=Peak_location[A$peak,1]
x1=Symbol_location[A$symbol,2]
y1=Peak_location[A$peak,2]
AA=(x==y)
BB=(abs(as.numeric(x1)-as.numeric(y1))<200000)
CC=abs(as.numeric(x1)-as.numeric(y1))
chr_idx=as.logical(AA*BB)
i=A$symbol[chr_idx]
p=A$peak[chr_idx]
d=CC[chr_idx]
x=rep(1,length(i))
beta_nonzero_idx=Matrix::sparseMatrix(p,i,x = x,dims=c(dim(Peak_location)[1],dim(Symbol_location)[1]))
beta_dis=Matrix::sparseMatrix(p,i,x = d,dims=c(dim(Peak_location)[1],dim(Symbol_location)[1]))
# filtet TGs without RE
TG_filter=TG_filter[Matrix::colSums(beta_nonzero_idx)>0,]
beta_dis=beta_dis[,Matrix::colSums(beta_nonzero_idx)>0]
beta_nonzero_idx=beta_nonzero_idx[,Matrix::colSums(beta_nonzero_idx)>0]
#tfidf
RE_filter@x=rep(1,length(RE_filter@x))
tf1=RE_filter/(matrix(1,ncol=1,nrow=dim(RE_filter)[1])%*%log(1+Matrix::colSums(RE_filter)))
idf=log(1+dim(RE_filter)[2])/(1+Matrix::rowSums(RE_filter))
O1=tf1*(as.matrix(idf)%*%matrix(1,nrow=1,ncol=dim(RE_filter)[2]));
O1[is.na(O1)]=0;
#impute RE TG matrix
library(Seurat)
O_s=CreateSeuratObject(O1)
O_s <- FindVariableFeatures(O_s, selection.method = "vst", nfeatures = floor(dim(RE_filter)[1]/4))
O_s <- ScaleData(O_s)
O_s <- RunPCA(O_s,npcs = 100)#!!! time consuming
score=O_s@reductions$pca@cell.embeddings
O_s=CreateSeuratObject(TG_filter)
O_s <- FindVariableFeatures(O_s, selection.method = "vst", nfeatures = floor(dim(TG_filter)[1]/4))
O_s <- ScaleData(O_s)
O_s <- RunPCA(O_s,npcs = 100)#!!! time consuming
score_TG=O_s@reductions$pca@cell.embeddings
TG_filter1=impute_dis(TG_filter,score_TG)
RE_filter1=impute_dis(O1,score)
RE_filter=O1
#initial
m=dim(RE_filter1)[2];
n=dim(TG_filter1)[2]#cell by gene
i=1:dim(RE_filter1)[1]
p=rep(1,dim(RE_filter1)[1])
x=rep(1,dim(RE_filter1)[1])
one_s=Matrix::sparseMatrix(i,p,x=x,dims=c(dim(RE_filter1)[1],1))
X=cbind(one_s,RE_filter1)
X1=as.matrix(X)
X=t(X1)%*%X1;#!!! time consuming
Y1=as.matrix(TG_filter1)
Y=t(Y1)%*%Y1
RE11=X1
b=exp(10*rowMeans(RE_filter1>0)); 
TG_mean=colMeans(TG_filter1);
TG_predict=as.matrix(b)%*%t(as.matrix(TG_mean))/mean(b);
beta = matrix(0,nrow=m,ncol=n);
TG_predict_meanZero=TG_predict-matrix(1,nrow=dim(TG_predict)[1],ncol=1)%*%t(as.matrix(TG_mean));
for (i in 1:n){
  X_temp=as.matrix(RE_filter1[,beta_nonzero_idx[,i]==1]);
  beta[beta_nonzero_idx[,i]==1,i] = pinv(t(X_temp)%*%X_temp)%*%t(X_temp)%*%TG_predict_meanZero[,i];
}
beta=rbind(t(TG_mean),beta)
beta_nonzero_idx=rbind(matrix(1,ncol=ncol(beta_nonzero_idx),nrow=1),beta_nonzero_idx)
beta_dis=rbind(matrix(0,ncol=ncol(beta_dis),nrow=1),beta_dis)
for (iter in 1:5){
  beta = Estimate_quadratic(X,Y,beta,beta_nonzero_idx,beta_dis,d0,lambda);
}
TG_predict1=matrix(0,nrow=dim(TG_predict)[1],ncol=ncol(TG_predict))

for (i in 1:n){
  TG_predict1[,i]=RE11[,beta_nonzero_idx[,i]==1]%*%beta[beta_nonzero_idx[,i]==1,i];
}
TG_predict3=TG_predict1*(b%*%matrix(1,nrow=1,ncol=dim(TG_predict)[2]))/mean(b)# %adjust by RE data
TG_predict3[TG_predict3<0]=0
colnames(TG_predict3)=rownames(TG_filter)
rownames(TG_predict3)=colnames(TG_filter)
beta1=beta[2:dim(beta)[1],]
# R=Matrix::t(TG_filter)%*%TG_filter;#TG-TG
R=t(TG_filter1)%*%(TG_filter1);#TG-TG
R=R*dim(TG_predict3)[1]/dim(TG_filter1)[1];
lambda=0.5;
peak_count=rowSums(O1>0);
f=order(peak_count);
n=floor(length(peak_count)/2);
R2=Matrix::t(O1[f[1:n],])%*%O1[f[1:n],]#use the half peaks highly specific
R2_rna=t(TG_filter1)%*%(TG_filter1);#cell by cell
scal=mean(R2_rna)/dim(TG_filter1)[1]*dim(O1)[2];
R2_rna=R2_rna/dim(TG_filter1)[1]*dim(O1)[2];
R2=R2/mean(R2)*scal;
d=diag(R2)/mean(diag(R2))*mean(diag(R2_rna));
R2=R2-diag(diag(R2))+diag(d);
X=TG_predict3;
R=as.matrix(R)
R2=as.matrix(R2)
for (i in 1:200){
  temp=(1+lambda)*X%*%t(X)%*%X
  temp[temp<10^(-20)]=10^(-20)
  X=X*(X%*%R+lambda*R2%*%X)/temp
}
rownames(beta1)=rownames(RE_filter)
return(list(X,beta1))
}
