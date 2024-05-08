clear, clc

n = 2;
u = zeros(n+1, n+1, n+1);
% [n1,n2,n3] = size(u);
% i = 2:n1-1;
% j = 2:n2-1;
% k = 2:n3-1;
% Au(i,j,k) = 6*u(i,j,k) - u(i-1,j,k) - u(i+1,j,k) - u(i,j-1,k) - u(i,j+1,k) ...
%     - u(i,j,k-1) - u(i,j,k+1);

e = ones(n+1,1);
D2 = spdiags([-e 2*e -e], -1:1, n+1, n+1);
E2 = kron(speye(n+1),D2) + kron(D2,speye(n+1));
F2 = kron(speye(n+1), E2) + kron(D2, kron(speye(n+1), speye(n+1)));
full(F2)