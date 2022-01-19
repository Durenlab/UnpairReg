function beta = Estimate_quadratic(X,Y,beta,beta_nonzero_idx,max_iter,enhancer_distance_for_gene,d0,lambda)
%%solving problem of Y=beta'*X*beta, subject to beta(~beta_nonzero_idx)=0;
n = size(Y,1);
%Xbeta=X*beta;
Xbeta=zeros(size(X,1),n);
for i=1:n
    Xbeta(:,i)=X(:,beta_nonzero_idx(:,i)==1)*beta(beta_nonzero_idx(:,i)==1,i);
end
for ii=1:max_iter
for i=1:n
TFAS=[0
exp(enhancer_distance_for_gene{i}/d0)];
W=TFAS*TFAS';
Temp=Xbeta(beta_nonzero_idx(:,i)~=0,setdiff([1:n],i))';
beta(beta_nonzero_idx(:,i)~=0,i)=pinv(Temp'*Temp+lambda*W)*Temp'*Y(setdiff([1:n],i),i);%pinv inverse
Xbeta(:,i)=X(:,beta_nonzero_idx(:,i)~=0)*beta(beta_nonzero_idx(:,i)~=0,i);
end
end
