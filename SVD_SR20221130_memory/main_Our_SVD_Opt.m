function [max_mem] = main_Our_SVD_Opt(a, c, r)

% % % 

n = size(a,1); 
d = sum(a,1)'; %  d: in-degree vector 
inv_d = spfun(@(x) 1./x, d);
q = a * spdiags(inv_d, 0, n, n)  ;   % q = col_norm(a)
max_mem = 0;
mem = whos;
max_mem = max(max_mem, sum([mem.bytes]));
clear d inv_d a ;
ide_n = speye(n);
[v, si, u] = svds(q, r);
mem = whos;
max_mem = max(max_mem, sum([mem.bytes]));
clear q n; 

% % % Improved Step:

g = v'* u * si; 
men = whos;
max_mem = max(max_mem, sum([men.bytes]));
clear v
ide_r = eye(r);
h = ide_r;
for i = 1:5
      h = c * g * h *  g' + ide_r; 
end
men = whos;
max_mem = max(max_mem, sum([men.bytes]));
clear g ide_r r
gamma = si * h * si;
men = whos;
max_mem = max(max_mem, sum([men.bytes]));
 clear si h

% % %

x = u * gamma;
men = whos;
max_mem = max(max_mem, sum([men.bytes]));
clear gamma
S_ap = (1-c)* ( ide_n + c * x * u');
% % %
men = whos;
max_mem = max(max_mem, sum([men.bytes]));
max_mem = max_mem/(1024^2);
fprintf('>>>>>>>>Algorithm 2(Our_SVD_Opt):  MaxMemory: %f MB\n ',  max_mem);

end

