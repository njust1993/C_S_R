function [max_mem] = main_Li_SVD( a, c, r)
% % % 
n = size(a,1); 
d = sum(a,1)'; %  d: in-degree vector 
inv_d = spfun(@(x) 1./x, d);
q = a * spdiags(inv_d, 0, n, n)  ;   % q = col_norm(a)
max_mem = 0;
mem = whos;
max_mem = max(max_mem, sum([mem.bytes]));
clear d inv_d a;   
ide = speye(n);
[v, si, u] = svds(q, r);     % w : svd-decomposition
mem = whos;
max_mem = max(max_mem, sum([mem.bytes]));
clear q r; 
% % %
K_u = kron(u, u);
mem = whos;
max_mem = max(max_mem, sum([mem.bytes]));

clear u;
inv_si= inv(si);
mem = whos;
max_mem = max(max_mem, sum([mem.bytes]));
clear si;

K_v = kron(v', v');
mem = whos;
max_mem = max(max_mem, sum([mem.bytes]));
clear v;
K_vu =  K_v * K_u;
lambda = inv(kron(inv_si, inv_si)-c * K_vu);  
mem = whos;
max_mem = max(max_mem, sum([mem.bytes]));
clear K_si  K_vu;
V_r = K_v * ide(:);
mem = whos;
max_mem = max(max_mem, sum([mem.bytes]));
clear K_v
% % %

P_r = K_u * lambda;
mem = whos;
max_mem = max(max_mem, sum([mem.bytes]));
clear K_u lambda;

S_ap = (1-c)*(ide(:)+c * P_r * V_r);
mem = whos;
max_mem = max(max_mem, sum([mem.bytes]));
max_mem = max_mem/(1024^2);
% % %
fprintf('>>>>>>>>>Algorithm 3( Li_SVD ):     MaxMemory: %f MB\n ',  max_mem);

end

