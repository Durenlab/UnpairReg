function [COEFF,SCORE,LATENT,EXPLAINED,COEFFss] = fastpca(data)
% -------------------------------------------------------------------------
% [COEFF,SCORE,LATENT,EXPLAINED] = fastpca(data) 
%
% Fast principal component analysis for very high dimensional data 
% (e.g. voxel-level analysis of neuroimaging data), implemented according
% to C. Bishop's book "Pattern Recognition and Machine Learning", p. 570.
% For high-dimensional data, fastpca.m is substantially faster than 
% MATLAB's in-build function pca.m.
% 
% According to MATLAB's PCA terminology, fastpca.m needs an input-matrix 
% with each row represents an observation (e.g. subject) and each column a 
% dimension (e.g. voxel). fastpca.m returns principal component (PC) 
% loadings COEFF, PC scores (SCORE), variances explained by the PCs in 
% absolute values (LATENT) and in percent (EXPLAINED). Additionally, 
% fastpca returns the PC loading of the small covariance matrix (COEFFs).
% 
% Decrease in computation time results from calculating the PCs from the 
% (smaller) covariance matrix of the transposed input-matrix "data" instead 
% of the large covariance matrix of the original input matrix which are 
% then use to project the observations to achieve the PCs of the large DxD 
% covariance matrix. 
% 
% By default, fastpca removes the mean of each observation.  In this first 
% implementation of fastpca, I skipped calculation of Hotelling痴 T-Squared 
% Statistic as I didn't need it so far.
%
% Example:
% In medical image analysis, there are often datasets with few to several 
% hundreds of observations (subjects) and hundreds of thousands dimensions 
% (voxels). As an example, I compare MATLABs PCA and fastpca using a random 
% matrix with 300 rows (e.g. subjects) and 500000 columns (e.g. voxels): 
%
%   data = rand(300,500000);
%
%   tic; [COEFF,SCORE,LATENT,~,EXPLAINED] = pca(data); toc
%   >> Elapsed time is 37.295108 seconds.
%
%   tic; [COEFF,SCORE,LATENT,EXPLAINED] = fastpca(data); toc
%   >> Elapsed time is 4.853614 seconds.
%
% Version 1.0 from 08/08/2019. 
% Implemented by Dominik Blum. 
% E-Mail: dominik.blum@med.uni-tuebingen.de
% Homepage: https://www.medizin.uni-tuebingen.de/de/das-klinikum/
%           mitarbeiter/profil/284?search=dominik%20Blum&mode=popup
% -------------------------------------------------------------------------
% 
[N,D] = size(data);
if N > D
    error('There are more observations than dimensions. fastpca.m is not suitable for this setting.')
end 
% centering observations
X = data - mean(data,1);
% eigenvalue decomposition of small covariance matrix
[COEFFs, LATENTs] = eig(1/(size(X,2)-1)*X*X');
[LATENTs,ind]     = sort(diag(LATENTs),'descend');
COEFFs            = COEFFs(:,ind);
% when observations are centered, the last PC will capture zero variance.
% If this is the case, this PC will be removed.
ind_neg           = find(LATENTs < 10^-10);
LATENTs(ind_neg)  = [];
COEFFs(:,ind_neg) = [];
% projection of scaled COEFFs onto observations to calculate PC loadings
COEFFss   = COEFFs*diag(1./sqrt((size(X,2)-1)*LATENTs));
COEFF     = X' * COEFFss; 
SCORE     = X * COEFF;
LATENT    = LATENTs.*((size(X,2)-1)/(size(X,1)-1));
EXPLAINED = 100 * LATENTs./sum(LATENTs);
