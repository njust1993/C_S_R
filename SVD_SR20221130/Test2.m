n = 1000;
a = spones(sprand(n, n, 0.9));
n = size(a,1); 
d = sum(a,1)'; %  d: in-degree vector 
inv_d = spfun(@(x) 1./x, d);
q = a * spdiags(inv_d, 0, n, n)  ;   % q = col_norm(a)
clear d inv_d a ;
 ide_n = speye(n);
r = 100;
[v, si, u] = svds(q, r);

g = v'* u * si; 
gt = g';
ide_r = eye(r);
x = ide_r;
for i = 1:25
      x = c * g * x * gt + ide_r;   
end

tic;
gamma = si * x * si;   %%% quckly
toc

tic;
xi = diag(si);
Gama = (xi*xi').* x;   %%% low
toc

