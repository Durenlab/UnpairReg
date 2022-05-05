# UnpairReg
UnpairReg deals with unpaired single cell multi-omics data, providing accurate estimation of cis-regulatory network and gene expression for cells in which chromatin accessibility is available.

**Run**

Input:

The input data could by an htf5 file including RNA-seq data and ATAC-seq data. We provide the code for PBMC 10x genomic data. We first download data:
```
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
```
Then we can download the code to the same folder 
```
git clone https://github.com/Durenlab/UnpairReg.git
```
Then, we can run UnpairReg in R:
```
source("h5tom.R")
source("impute_dis.R")
source("UnpairReg.R")
source("Estimate_quadratic.R")
result=h5tom("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
predict=UnpairReg(result)
TG_prediction=predict[[1]]
cis_regulate=predict[[2]]
```
Output:

The output including TG_prediction and cis_regulate, which represent the predicted gene expression and the cis-regulatory network.


**Requirements**

R packages: rhdf5, stringr, tidyr, Matrix, pracma, Seurat
