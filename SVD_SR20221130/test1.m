n = 1000;
a = spones(sprand(n, n, 0.9));
c = 0.6;
r = 100;
n = size(a,1); 
d = sum(a,1)'; %  d: in-degree vector 
inv_d = spfun(@(x) 1./x, d);
q = a * spdiags(inv_d, 0, n, n)  ;   % q = col_norm(a)
[v, si, u] = svds(q, r);
% inv_si = spdiags(1./diag(si),0, r, r);
inv_si = inv(si);   % The dimension of si is less than 100 and the method is faster 
theta = v' * u;
lambda = inv(kron(inv_si, inv_si) - c * kron(theta, theta ));
ide = eye(r);
gamma = reshape (lambda * ide(:), [r,r]);

%%%%%%%%%%%%
a1 = v'* u * si;  
a2 = a1';
x1 = ide;

ErrorBound = 1.0e-4; % % % or ErrorBound = 1.0e-5;
kmax =  round(log( ErrorBound)/log(c));
% 老师，如果设置 c = 0.8时；若使误差保持在 1.0e-6 左右，迭代次数至少在 kmax在 52 - 62 左右；  若误差保持在 1.0e-4  左右，迭代次数 kmax在 31 - 41 左右； 
% 老师，如果设置 c = 0.6 时；若使误差保持在 1.0e-6 左右，迭代次数至少在 kmax在23 - 27 左右；  若误差保持在 1.0e-4  左右，迭代次数 kmax在 14-18 左右； 

for k1 = 1: kmax 
x1 = c* a1 * x1 * a2 + ide;
end

Ga = si * x1 * si;
sum(sum(gamma-Ga))