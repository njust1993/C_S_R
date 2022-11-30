function [max_mem] = I_Sim(a, c, kmax)   % s = cw'sw + (1-c)I
%%%%%%  computate the peak amount of memory used by the algorithm to run.


n = size(a,1); 
d = sum(a,1)'; %  d: in-degree vector 
inv_d = spfun(@(x) 1./x, d);
w = a * spdiags(inv_d, 0, n, n)  ;   % q = col_norm(a)
max_mem = 0;
mem = whos;
max_mem = max(max_mem, sum([mem.bytes]));
clear  inv_d a d;


ide = speye(n, n);
s = ide;
mem = whos;
max_mem = max(max_mem, sum([mem.bytes]));
clear n

%%%
for i = 1: kmax
    s = c * w' * s * w + (1-c) * ide;
    if(eq(i, 1))
        men = whos;
        max_mem = max(max_mem, sum([men.bytes]));
    end
    
end
max_mem = max_mem/(1024^2);
fprintf('>>>>>>>>>>Algorithm 1: Iterative method of SimRank(I_Sim) (baseline)                          MaxMemory: %f MB\n ',  max_mem); 
fprintf('\n');
end

