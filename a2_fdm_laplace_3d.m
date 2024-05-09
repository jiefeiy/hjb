clear, clc
% solve 3D Poisson's equation by finite difference method
%
% -\Delta u = f,   in [0,\pi]^3,
% with periodic boundary condition, 
% where f(x,y) = 8*pi^2*sin(2*pi*X).*cos(4*pi*Y).*cos(2*pi*Z);   % source term
%
% first generate a sparse 3D negative Laplacian matrix
% define component matrices
n = 30;
e = ones(n+1,1);
D1x = spdiags([-e 2*e -e], [-1 0 1], n+1, n+1);
D1y = D1x; D1z = D1x;

% set period boundary condition
for i = 1:3
    eval(['D1' char(119 + i) '(1,' num2str(n+1) ') = D1'...
                char(119 + i) '(1,' num2str(n+1) ') -1;']);
    eval(['D1' char(119 + i) '(' num2str(n+1) ',1) = D1'...
                char(119 + i) '(' num2str(n+1) ',1) -1;']);
end

% generate matrix by Kronecker tensor product
A = kron(speye(n+1), kron(speye(n+1), D1x)) ...
    + kron(speye(n+1), kron(D1y, speye(n+1))) ...
    + kron(kron(D1z, speye(n+1)), speye(n+1));

% solving linear system
h = 1/n;
[X, Y, Z] = meshgrid(0:h:1, 0:h:1, 0:h:1);
f = h^2*8*pi^2*sin(2*pi*X).*cos(4*pi*Y).*cos(2*pi*Z);       % source term * h^2
b = reshape(f,[],1); 
u_vec = A \ b;
u = reshape(u_vec, [n+1, n+1, n+1]);

surf(X(:,:,1), Y(:,:,1), u(:,:,1));
