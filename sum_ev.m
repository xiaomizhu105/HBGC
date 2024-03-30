function [sum_k,sum_k_1] = sum_ev(P,k)
[n,m] = size(P);
P = sparse(P);

p1 = sum(P,2);
Dn = spdiags(1./sqrt(p1),0,n,n);

p2 = sum(P,1);
Dm = spdiags(1./sqrt(p2'),0,m,m);

S = sparse(n+m,n+m); 
S(1:n,n+1:end)=P; 
S(n+1:end,1:n)=P';

Ds = sparse(n+m,n+m); 
Ds(1:n,1:n) = Dn;
Ds(n+1:end,n+1:end) = Dm;

Ls = eye(n+m,n+m) - Ds*S*Ds; % n+m x n+m
Ls = sparse(Ls);
ev = svds(Ls,k+1,'smallest');
sum_k = sum(ev(2:end));
sum_k_1 = sum(ev);
end

