function [Total_Time] = main_Li_SVD( a, c, r, s, Time_Pro_SVD)
     % % % 
% Pro_SVD = tic;   
n = size(a,1); 
d = sum(a,1)'; %  d: in-degree vector 
inv_d = spfun(@(x) 1./x, d);
q = a * spdiags(inv_d, 0, n, n)  ;   % q = col_norm(a)
clear d inv_d a;   
ide = speye(n);
[v, si, u] = svds(q, r);     % w : svd-decomposition
% Time_Pro_SVD = toc(Pro_SVD);
clear q; 

 % % %
InterPro= tic;
K_u = kron(u, u);
clear u;
inv_si= inv(si);
clear q si;
K_v = kron(v', v');
clear v;
K_vu =  K_v * K_u;
lambda = inv(kron(inv_si, inv_si)-c * K_vu);     
clear K_si  K_vu;
V_r = K_v * ide(:);
Time_InterPro = toc( InterPro); 

% % %
Com_Sim = tic;
P_r = K_u * lambda;
clear K_u lambda;
S_ap = (1-c)*(ide(:)+c * P_r * V_r);  
Time_Com_Sim = toc( Com_Sim );
% % %

AvgDiff = sum(abs(S_ap - s(:)))/(n^2);       %%%  Average Difference(EDBT)
Total_Time = Time_Pro_SVD + Time_InterPro + Time_Com_Sim;
fprintf('>>>>>>>>>Algorithm 3( Li_SVD ):    Total_Time = %f s; Time_Pro_SVD = %f; Time_InterPro = %f;  Time_Com_Sim = %f;  AvgDiff = %f;\n', Total_Time, Time_Pro_SVD, Time_InterPro, Time_Com_Sim, AvgDiff);

end

