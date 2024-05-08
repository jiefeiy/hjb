clear, clc
% solve 2D Poisson's equation by finite difference method
%
% -\Delta u = f,   in [0,1]\times[0,1],
%    u(x,1) = 2,   on 0\le x \le 1,
% where f(x,y) = 8*pi^2*sin(2*pi*X).*cos(2*pi*Y);   % source term
% 
% Reference: 
% [1] https://www.math.uci.edu/~chenlong/226/FDM.pdf
% [2] https://www.math.uci.edu/~chenlong/226/FDMcode.pdf

%% 2d finite difference
upbound = 2;
n = 40;
h = 1/n;
[X, Y] = meshgrid(0:h:1,0:h:1);
f = h^2*8*pi^2*sin(2*pi*X).*cos(2*pi*Y);       % source term * h^2

%% solve only equations at interior nodes
e = ones(n-1,1);
T = spdiags([-e 2*e -e], -1:1, n-1, n-1);      % matrix for pts in the middle
A = kron(speye(n-1),T) + kron(T,speye(n-1));

b = reshape(f(2:n,2:n)',[],1); 
b((n-1)*(n-2)+1:end) = b((n-1)*(n-2)+1:end) + upbound;  % boundary condition, upper boundary = 2
u_vec = A\b; 

%% display solution
u = zeros(n+1);
u(n+1, :) = upbound;                           % boundary condition
u(2:n, 2:n) = reshape(u_vec, n-1, n-1)';
% contourf(X, Y, u, 10);
surf(X,Y,u);
xlabel('x'); ylabel('y'); zlabel('u');
