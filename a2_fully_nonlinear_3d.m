clear, clc
% numerical tests:
% [1] W Guo, J Zhang, J Zhuo. (2015) A monotone scheme for high-dimensional
% fully nonlinear pdes.
%
% test with example 1 in Section 6 Numerical examples. 
%
% 

%% parameter setting 
% sigma_min = 1;
% sigma_max = sqrt(2);
% T = 0.5;
% x0 = [5 6 7]';

%% 2d finite difference
n = 4;
e = ones(n,1);
T = spdiags([-e 2*e -e], -1:1, n, n);


