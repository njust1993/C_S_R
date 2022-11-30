function [s] = I_Sim(a, c, kmax)   % s = cw'sw + (1-c)I
%%%%%%  computate the peak amount of memory used by the algorithm to run.

time_Iter = tic;
n = size(a,1); 
d = sum(a,1)'; %  d: in-degree vector 
inv_d = spfun(@(x) 1./x, d);
w = a * spdiags(inv_d, 0, n, n)  ;   % q = col_norm(a)
clear  inv_d a d;
wt = w';

ide = speye(n, n);
s = ide;
clear n

%%%
for i = 1: kmax
    s = c * wt * s * w + (1-c) * ide;    
end

Temp_Iter = toc(time_Iter);
fprintf('>>>>>>>>>>Algorithm 1: Iterative method of SimRank(I_Sim) (baseline)                          : Time: %f s\n ',Temp_Iter); 
fprintf('\n\n');
end

