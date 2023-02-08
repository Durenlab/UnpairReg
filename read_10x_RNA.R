read_10x_RNA=function(dir,genome){ 
library(Seurat)

# Access the count matrix
E <- Read10X(dir)
#gene location Human
if (genome=='hg19'){
    gene_loc=read.table('hg19.gene_anno.GeneBreak.txt',sep='\t',header=FALSE)
    colnames(gene_loc)=c('Gene','EnsID','Chromosome','Start','End','band','strand')
}else if(genome=='hg38'){
gene_loc=read.table('hg38_gene_annotation.txt',sep='\t',header=TRUE)
    }else{print('genome is not accepted')}
Symbol=rownames(E)
loc=match(gene_loc$Gene,Symbol)
gene_loc=gene_loc[!is.na(loc),]
tss=gene_loc$Start
tss[gene_loc$strand==-1]=gene_loc$End[gene_loc$strand==-1]
E=E[loc[!is.na(loc)],]
Symbol_location=cbind(gene_loc$Chromosome,tss)
return(list(E,Symbol_location)) 
}
