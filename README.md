# UnpairReg
UnpairReg deals with unpaired single cell multi-omics data, providing accurate estimation of cis-regulatory network and gene expression for cells in which chromatin accessibility is available.

# Run
Input:
1. Single cell RNA-seq count matrix, 'E'. Rows respresent cell and columns respresent gene;
2. Sinle cell ATAC-seq count matrix, 'O'. Rows respresent cell and columns respresent regulatory element (RE); 
3. Gene symbol of RNA-seq data, 'Symbol', which is a matlab cell vector;
4. Gene transcriptional start site (TSS), a matrix with 2 columns, 'Symbol_location'. The first is the number of chromosome ('X' as 23 and 'Y' as 24); the second is the location on chromosome.
5.  Location of  RE, 'Peak_location', whose format is same as 'Symbol_location'.
All inputs above are embedded into 'Input.mat'.

Run:
matlab -nodisplay -nosplash -nodesktop -r "code; exit"
