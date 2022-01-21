function run(k,lambda,range,d0,outdir)
load('Input.mat','E','O','Symbol_location','Symbol','Peak_location')
RE_filter=full(log(O+1));
gene_loc=Symbol_location;
enhancer_loc=Peak_location;
enhancer_index_for_gene = cell(size(gene_loc,1),1);%save the enhancer for the ith gene within 1M
enhancer_distance_for_gene = cell(size(gene_loc,1),1);
index_all=1:size(enhancer_loc,1);
for i = 1:size(gene_loc,1)
    temp_chr = gene_loc(i,1);
    temp_loc = gene_loc(i,2);
    index_chr = (enhancer_loc(:,1) == temp_chr);
    index1 = index_all(index_chr);
    enhancer_loc_temp_chr = enhancer_loc(index_chr,2);
    index_loc = (abs(enhancer_loc_temp_chr-temp_loc)<range);
    enhancer_index_for_gene{i}=index1(index_loc);
    enhancer_distance_for_gene{i}=abs(enhancer_loc_temp_chr(index_loc)-temp_loc);
end
beta_nonzero_idx=zeros(size(O,2),size(E,2));
for i=1:length(enhancer_index_for_gene)
    beta_nonzero_idx([1,enhancer_index_for_gene{i}+1],i)=1;
end
O=1*(O'>0);
tf1=O./(ones(size(O,1),1)*log(1+sum(O)));
idf=log(1+size(O,2)./(1+sum(O>0,2)));
O1=tf1.*(idf*ones(1,size(O,2)));
O1(isnan(O1))=0;
RE=full(O1');
[coeff0,score,latent] = fastpca(RE);
clear coeff0 latent;
dim100_RE=score(:,1:100);
TG1=full(log(E+1)); %observed TG
[coeff0,score,latent] = fastpca(TG1);
clear coeff0 latent;
dim100_TG=score(:,1:100);
clear KNN1 KNN2;
clear O tf1 idf O1;
%%%%%%%%%%%%%%%%imputation
TG_filter=full(TG1);
beta_nonzero_idx=full(beta_nonzero_idx);
[RE_filter1,KNN1]=immupute_dis(RE_filter,dim100_RE,k);
[TG_filter1,KNN2]=immupute_dis(TG_filter,dim100_TG,k);
clear KNN1 KNN2;
m=size(RE_filter1,2);
n=size(TG_filter1,2);
X=[ones(size(RE_filter1,1),1),RE_filter1]'*[ones(size(RE_filter1,1),1),RE_filter1];
Y=TG_filter1'*TG_filter1;
beta_nonzero_idx=logical(beta_nonzero_idx);
%%% initial of beta
RE11=[ones(size(RE_filter1,1),1),RE_filter1];
b=exp(10*mean(RE_filter1'>0)'); %%seq_depth
TG_mean=mean(TG_filter1);
TG_predict=b.*TG_mean/mean(b);
beta = zeros(m,n);
TG_predict_meanZero=TG_predict-TG_mean;
for i=1:n
	X_temp=full(RE_filter1(:,beta_nonzero_idx(2:end,i)==1));
	beta(beta_nonzero_idx(2:end,i)==1,i) = pinv(X_temp'*X_temp)*X_temp'*TG_predict_meanZero(:,i); % gain the beta by linear regression for each gene with the nearby enhancer
end
beta=[TG_mean;beta];
%%%%%%%%%%%%%%%%%%%%Estimate beta
for iter=1:2
beta = Estimate_quadratic(X,Y,beta,beta_nonzero_idx,1,enhancer_distance_for_gene,d0,lambda); % solve beta 32522.317481(1M) 9hours;200k(2h)
for i=1:n
    TG_predict(:,i)=RE11(:,beta_nonzero_idx(:,i)==1)*beta(beta_nonzero_idx(:,i)==1,i);
end
TG_predict=max(TG_predict,0);
b=exp(10*mean(RE_filter1'>0)');
TG_predict3=TG_predict.*b/mean(b); %adjust by RE data
end
%%output
save('Output.mat','beta','TG_predict3','Symbol','-v7.3');
