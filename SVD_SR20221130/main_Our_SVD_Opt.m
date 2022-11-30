function [Time_Pro_SVD] = main_Our_SVD_Opt(a, c, r, s)

% % % 
Pro_SVD = tic; 
n = size(a,1); 
d = sum(a,1)'; %  d: in-degree vector 
inv_d = spfun(@(x) 1./x, d);
q = a * spdiags(inv_d, 0, n, n)  ;   % q = col_norm(a)
clear d inv_d a ;
 ide_n = speye(n);
[v, si, u] = svds(q, r);
Time_Pro_SVD = toc(Pro_SVD);
clear q ; 

% % % Improved Step:
InterPro= tic;
g = v'* u * si; 
gt = g';
ide_r = eye(r);
x = ide_r;
for i = 1:5
      x = c * g * x * gt + ide_r;   
end
gamma = si * x * si;
Time_InterPro = toc( InterPro); 

% % %
Com_Sim = tic;
x = u * gamma;
clear gamma
S_ap = (1-c)* ( ide_n + c * x * u');
% % %
Time_Com_Sim = toc( Com_Sim );

AvgDiff = sum(sum(abs(S_ap-s)))/(n^2);       %%% Average Difference
Total_Time = Time_Pro_SVD + Time_InterPro + Time_Com_Sim;

fprintf('>>>>>>>>Algorithm 2(Our_SVD_Opt):  Total_Time = %f s; Time_Pro_SVD = %f; Time_InterPro = %f;  Time_Com_Sim = %f;  AvgDiff = %f;\n', Total_Time, Time_Pro_SVD, Time_InterPro, Time_Com_Sim, AvgDiff);

end

