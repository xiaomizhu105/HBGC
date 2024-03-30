function [Fn,Fm,Dn,Dm] = max_Tr_FLsF(P,k)
[n,m] = size(P);
P = sparse(P);
p1 = sum(P,2);
Dn = spdiags(1./sqrt(p1),0,n,n);
p2 = sum(P,1);
Dm = spdiags(1./sqrt(p2'),0,m,m);
X = Dn*P*Dm; % n x m

[U,S,V] = svds(X,k);
Fn = sqrt(2)/2*U; 
Fm = sqrt(2)/2*V;
Dn = full(Dn);
Dm = full(Dm);
end

