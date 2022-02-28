DataInput('pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5');
k=100;
lambda=10^7;
range=200000;
d0=10000;
outdir='./';
run(k,lambda,range,d0,outdir);
