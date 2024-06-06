clear, clc, close all

% consider 1D Darcy flow problem 
%       - d/dx ( a(x,theta) d/dx(u(x,theta)) ) = f(x),  for x \in (0,1)
%       u(0,theta) = 0,   u(1,theta) = 0
% where 
%       f(x) = 1000 if x \in [0, 0.5]; 
%       f(x) = 2000 if x \in [0.5, 1].
% 
% The equation describes the pressure field u(x,theta) in a porous medium 
% with a random permeability field a(x, theta). 
% Here, x is the location and theta is the parameter to model randomness.
% a(x, theta) is modeled by a log-normal random field.
%
% Example: 
%   The quantity of interest is E[u(1/3, theta)], where the expectation is
%   taken over the parameter theta (a Nl-dimensional Gaussian vector). 


%% parameters
N = 512;
Nl = 32;
d = 2;
tau = 3;
f = @(x) 1000*(x<=0.5) + 2000*(x>0.5);         % source of fluid 
% f = @(x) ones(size(x));

%% KL expansion 
phi = @(l,x) sqrt(2)*cos(pi*l.*x);
lambda = @(l) (pi^2*l.^2 + tau.^2).^(-d);
sqrt_lambda_all = sqrt(lambda(1:Nl))';

theta_all = randn(Nl, 1);                      % randomness comes from here
a = @(x) exp( sum(theta_all .* sqrt_lambda_all .* phi((1:Nl)', x), 1) );
xgrid = linspace(0, 1, N+1);

%%% plot one realization of permeability field a(x,theta)
plot(xgrid, a(xgrid));
xlabel('x'); ylabel('permeability');

%% solve the pde for one realization of a(x,theta)
[~, u] = solve_darcy_1d(1/3, theta_all, f, N, Nl, d, tau);
figure(2);
plot(xgrid, u);
xlabel('x'); ylabel('u(x)');


%% solution u with respect to theta
% M = 100;                      % M collocation points of theta
% theta1 = linspace(-3,3,M);
% ux = zeros(M,1);
% 
% for m = 1:M
%     theta_all = [theta1(m); zeros(Nl-1, 1)];
%     ux(m) = solve_darcy_1d(1/3, theta_all, f, N, Nl, d, tau);
% end
% plot(theta1, ux);



%% dependencies
function [u_x, u] = solve_darcy_1d(x, theta, f, N, Nl, d, tau)
% solve the pde for one realization of a(x,theta)
% using finite difference scheme
phi = @(l,x) sqrt(2)*cos(pi*l.*x);
lambda = @(l) (pi^2*l.^2 + tau.^2).^(-d);
sqrt_lambda_all = sqrt(lambda(1:Nl))';

a = @(x) exp( sum(theta .* sqrt_lambda_all .* phi((1:Nl)', x), 1) );
xgrid = linspace(0, 1, N+1);

dx = 1/N;
agrid = a(xgrid);
A = diag( agrid(1:N-1) + 2*agrid(2:N) + agrid(3:N+1) );
A = A + diag( -(agrid(2:N-1) + agrid(3:N)) , -1) ...
    + diag( -(agrid(2:N-1) + agrid(3:N)) , 1);
fgrid = f(xgrid(2:N))' * (2*dx*dx);
u = A \ fgrid;
u = [0; u; 0];
u_x = interp1(xgrid, u, x, 'linear');
end










