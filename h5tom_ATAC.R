h5tom_ATAC=function(h5file){ 
library(rhdf5)
library(stringr)
data0=h5read(h5file,'/matrix/data')
indices=h5read(h5file,'/matrix/indices')
indptr=h5read(h5file,'/matrix/indptr')
shape=h5read(h5file,'/matrix/shape')
name=h5read(h5file,'/matrix/features/name')
barcode=h5read(h5file,'/matrix/barcodes')
i = as.integer(indices)
p = as.integer(indptr)
x = as.numeric(data0)
rp_matrix1=Matrix::sparseMatrix(c(1,3:8), c(2,9,6:10), x = 7 * (1:7),dims=shape)
rp_matrix1@i=i
rp_matrix1@p=p
rp_matrix1@x=x
colnames(rp_matrix1)=barcode
rownames(rp_matrix1)=name
name_s=str_split(name, ':',simplify = TRUE)
npeak=dim(name_s)[1]
O=rp_matrix1[1:dim(name_s)[1],]
#region=str_split(interval,c(':'),simplify = TRUE)
loc=str_split(name_s[,2],'-',simplify=TRUE)
chro=paste0('chr',c(1:22,'X','Y'))
ATAC_check=match(name_s[,1],chro)
A=1:npeak
id2=A[!is.na(ATAC_check)]
O=O[id2,]
Peak_location=cbind(name_s[1:dim(name_s)[1],1],loc[1:dim(name_s)[1],1])[id2,]
return(list(O,Peak_location))
}
