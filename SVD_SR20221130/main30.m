clear all
clc
%%%%  example 2 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 'email-Eu-core';              %% nodes: 1005,   edges: 25571                        % % N_Sim: max(r)=32      SVD_SR: max(r)=213 
% %  'ego-Facebook';               %% nodes: 4039,   edges: 88234                       % % N_Sim: max(r)=8       SVD_SR: max(r)=213 
% %  'ca-GrQc';                    %% nodes: 5,242,  edges:28980                        % % N_Sim: max(r)=5       SVD_SR: max(r)=212 
% %  'as-735';                     %% nodes: 7716,   edges:26467                        % % N_Sim: max(r)=5       SVD_SR: max(r)=210 
% %  ds = 'p2p-Gnutella08';        %% nodes: 6301    edges:20777                        % % N_Sim: max(r)=4       SVD_SR: max(r)=210 
% %  'wiki-Vote';                  %% nodes: 8297,   edges:103689                       % % N_Sim: max(r)=3       SVD_SR: max(r)=210 
% %  'p2p-Gnutella06';             %% nodes: 8717,    edges: 31525                      % % N_Sim: max(r)=3       SVD_SR: max(r)=210 
% %  'ca-HepPh';                   %% nodes: 12008,  edges:237010                       % % N_Sim: max(r)=2       SVD_SR: max(r)=210 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% 
 fpath = 'D:\Scientific research\dataset\';
%  fpath = 'E:\MATLAB\dataset\';
%1 % ds = 'email-Eu-core';                    %% nodes: 1005,   edges: 25571  
%1 % ds = 'Roget';                             %% nodes: 1022,    edges: 5075         rank =984
% %*  ds = 'ego-Facebook';                     %% nodes: 4039,   edges: 88234        r = 10 avgdiff (*)
%1 %  ds = 'power';                     %% nodes: 4941,   edges: 13188 
% % ds = 'USpowerGrid';                     %% nodes: 4941,   edges: 13188  % meixiazai 
% %* ds = 'delaunay_n13';                     %% nodes: 4941,   edges: 13188   % 内存不足
% %*  ds = 'ca-GrQc';                          %% nodes: 5242,      edges:28980      r = 10 Time(*) 
% %*  ds = 'Erdos972';                         %% nodes: 5488,      edges:14170 * *   r = 10 Time(*)    
% %*  ds = 'p2p-Gnutella08';                   %% nodes: 6301       edges:20777 *
% %*  ds = 'delaunay_n11';                     %% nodes: 2048       edges:12032 *
%1 %  ds = 'as-735';                           %% nodes: 7716,      edges:26467    
%1 %  ds = 'wiki-Vote';                        %% nodes: 8297,      edges: 103689 
% %  ds = 'p2p-Gnutella06';                   %% nodes: 8717,      edges: 31525
% %  ds = 'ca-HepTh';                         %% nodes: 9877,      edges:51971 *
% %* ds = 'ca-HepPh';                         %% nodes: 12008,     edges:237010   *
%1 %  ds = 'ca-AstroPh';                       %% nodes: 18772,     edges:396,160
%1 %  ds = 'p2p-Gnutella25';                   %% nodes: 22687,     edges: 54705    
%1 %  ds = 'cit-HepPh';                        %% nodes: 34546,     edges: 421578 
% %*  ds = 'email-Enron';                      %% nodes: 36692,     edges: 367662        %内存不足
% %*  ds = 'web-NotreDame';                    %% nodes: 325729,     edges: 1497134   %内存不足


  a = loaddata(fpath, ds);
      
 
         
%%%% example 3 
%%%%  stochastic matrix
% n = 20;
% a = spones(sprand(n, n, 0.9));

c = 0.6;
 %%%%%% Algorithm 1: Iterative method (baseline):   %%  s = cw'sw + (1-c)I
 ErrorBound = 1.0e-4; % % % or ErrorBound = 1.0e-5;
 kmax =  round(log( ErrorBound)/log(c));
 [s] = I_Sim(a, c, kmax);
 
 
 r = [5 10 50 100 200];      % r can change, r-svd decomposition                          
%  r = [5 10 50];                           
 parameter = 3;      

 for t = 1:size(r,2)
     
%%%%%% Algorithm 2: Our algorithm:  
fprintf('>>>>>>>>> r = %d;  Similarity Value : 1) time 2)avgdiff \n ', r(t));
fprintf('> \n');
[Time_Pro_SVD] = main_Our_SVD_Opt(a, c, r(t), s);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Algorithm 3: Li et al's algorithm (N_Sim)
%      if r(t) <= parameter
%             [Total_Time] = main_Li_SVD( a, c, r(t), s, Time_Pro_SVD); 
%      else
%             fprintf('> \n');
%             continue;
% 
%      end
%  fprintf('> \n');
 end