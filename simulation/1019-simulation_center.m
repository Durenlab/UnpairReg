cd /data/duren_lab/Kaya/Unpair/simulation/center/
k=100;
n=1000;
m=5000;
c=10000;
SNR=50;
K=floor(m*0.1);
centers=rand(10,m);
RE=zeros(c,m);
for i =1:10
Z = randn(1000, m);
RE((1000*(i-1)+1):1000*i,:)=centers(i,:)+1/10*Z;
end
RE=RE-min(min(RE));

scATAC_idx = randperm(c,floor(0.5*c)); % floor(0.5*c) randome unique integers from 1:c
scRNA_idx =  setdiff([1:c],scATAC_idx); %
beta_true = sprandn(m,n,10/m);
i=0;
j=0;
out=zeros(20,20);
for drop_out_rate_RE = linspace(0,0.95,20)
i=i+1;
j=0;
	for drop_out_rate_TG = linspace(0,0.95,20)
	j=j+1;
	tic 
[TG,TG1,RE1,beta_nonzero_idx]=data_simulation(n,m,c,RE,beta_true,SNR,drop_out_rate_RE,scATAC_idx,scRNA_idx,drop_out_rate_TG);
[coeff0,score,latent] = pca(RE1);
clear coeff0 latent;
dim100_RE=score(:,1:100);
[coeff0,score,latent] = pca(TG1);
clear coeff0 latent;
dim100_TG=score(:,1:100);
[TG_filter1,KNN2]=immupute_dis(TG1,dim100_TG(:,1:20),k);
[RE_filter1,KNN2]=immupute_dis(RE1,dim100_RE(:,1:20),k);
% 1/SNR is the scale of error matrix
X2 = RE_filter1'*RE_filter1; % 
Y2 = TG_filter1'*TG_filter1; 
%%% estimation
id = find(beta_nonzero_idx>0);
b=exp(10*mean(RE_filter1'>0)'); %%seq_depth
TG_mean=mean(TG_filter1);
TG_predict=b*TG_mean/mean(b);
beta = zeros(m,n);
for i=1:n
	X_temp=full(RE_filter1(:,beta_nonzero_idx(:,i)==1));
	beta(beta_nonzero_idx(:,i)==1,i) = pinv(X_temp'*X_temp)*X_temp'*TG_predict(:,i); % gain the beta by linear regression for each gene with the nearby enhancer
end
beta = Estimate_quadratic(X2,Y2,beta,beta_nonzero_idx,10,beta_true,1000); % solve beta
% evaluation
id = find(beta_nonzero_idx ~= 0);
full_beta_true=full(beta_true(id));
full_beta=full(beta(id));
toc
out(i,j)=corr(full_beta_true,full_beta);
end
end

save('data_normal.mat','scATAC_idx','scRNA_idx','beta_true','RE','centers')
