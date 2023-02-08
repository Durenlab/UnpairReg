# UnpairReg
---
UnpairReg deals with unpaired single cell multi-omics data, providing accurate estimation of gene expression for cells in which chromatin accessibility is available as well as cis-regulatory network.

The input is scRNA-seq and scATAC-seq data from the same tissue/context but from different cells. 

## Run

For the first step, we can download the code by
```
git clone https://github.com/Durenlab/UnpairReg.git
```
### Input:

The input data could be an htf5 file for ATAC-seq data, as well as RNA-seq data including barcodes.tsv, genes.tsv, and matrix.mtx. We provide the code for PBMC 10x genomic data. We first download data:
```sh
wget https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_pbmc_5k/atac_v1_pbmc_5k_filtered_peak_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-exp/1.1.0/frozen_bmmc_healthy_donor1/frozen_bmmc_healthy_donor1_filtered_gene_bc_matrices.tar.gz
gzip -d frozen_bmmc_healthy_donor1_filtered_gene_bc_matrices.tar.gz
```
Then, we can run UnpairReg in R:
```R
source('h5tom_ATAC.R')
source('read_10x_RNA.R')
h5file_ATAC='atac_v1_pbmc_5k_filtered_peak_bc_matrix.h5'
genome='hg19'
dir='filtered_matrices_mex/hg19/'#the dir incuding RNA-seq data
ATAC=h5tom_ATAC(h5file_ATAC)
RNA=read_10x_RNA(dir,genome)
RNA=read_10x_RNA(dir,genome)
result=list(RNA[[1]],ATAC[[1]],RNA[[2]],ATAC[[2]])
source("impute_dis.R")
source("UnpairReg.R")
source("Estimate_quadratic.R")
predict=UnpairReg(result,d0=10000,lambda=10^7)
TG_prediction=predict[[1]]
cis_regulate=predict[[2]]
```
### Input a multimodal data
```sh
cd UnpairReg
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
```

Then, we can run UnpairReg in R:
```R
source("h5tom.R")
source("impute_dis.R")
source("UnpairReg.R")
source("Estimate_quadratic.R")
result=h5tom("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
predict=UnpairReg(result,d0=10000,lambda=10^7)
TG_prediction=predict[[1]]
cis_regulate=predict[[2]]
```
### Output:

The output includes TG_prediction and cis_regulate, which represent the predicted gene expression and the cis-regulatory network.


### Requirements

R packages: rhdf5, stringr, tidyr, Matrix, pracma, Seurat. 
## Reference
Yuan, Qiuyue, and Zhana Duren. "Integration of single-cell multi-omics data by regression analysis on unpaired observations." Genome biology 23.1 (2022): 1-19. https://link.springer.com/article/10.1186/s13059-022-02726-7 
