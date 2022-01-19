function [X,KNN]=immupute_dis(X,H,k)
sim=squareform(pdist(H));
sim=max(max(sim))-sim;
%sim=sim.*(eye(size(sim,1))--0);
KK = full(sum(sim));
twom = sum(KK);
sim_norm=sim - KK'*KK/twom;
[d f]=sort(sim_norm,'descend');
KNN_P=[];
for i=1:size(sim_norm,2)
KNN_P=[KNN_P;[i*ones(k,1) f(1:k,i)]];
end
KNN=sparse(KNN_P(:,2),KNN_P(:,1),(1/k)*ones(size(KNN_P,1),1),size(sim_norm,2),size(sim_norm,2));
X=KNN'*X;
