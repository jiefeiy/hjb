clear, clc
% implement LQG control problem for d=1
% naive implement using finite difference + forward Euler
T = 1;
lambda = 2;
g = @(x) log((1 + x.^2)/2);  % terminal condition 

% xgrid = linspace(-100,100,100);
% plot(xgrid, g(xgrid), '.-');

% decide boundary condition 
xmin = -10;
xmax = -xmin;

% discretization
n_time = 50;
n_x = 50;
dx = (xmax - xmin)/n_x;
dt = T/n_time;
xgrid = (xmin:dx:xmax)';

u = zeros(n_x+1, n_time+1);
u(:, end) = g(xgrid);
for n = n_time:-1:1
    for j = 2:n_x
        u(j,n) = u(j,n+1) + dt/(dx^2)*(u(j-1,n+1) + u(j+1, n+1) - 2*u(j,n+1)) ...
            - lambda/4*dt/(dx)^2 * (u(j+1,n+1) - u(j-1,n+1)).^2;
    end
    u(1,n) = u(1,n+1);
    u(n_x+1, n) = u(n_x+1,n+1);
end

[X,Time] = ndgrid(xgrid, 0:dt:T);
surf(X, Time, u);

%% exact solution approximation 
% u(0,x) = -1/lambda * log( E[ ((1 + (x + sqrt(2)*sqrt(T)*randn)^2)/2)^(-lambda) ] )
rng("default")
M = 40000;
u0_exact = zeros(size(xgrid));
for j = 1:n_x+1
    u0_exact(j) = -1/lambda * log( mean( ((1 + (xgrid(j) + sqrt(2)*sqrt(T)*randn(M,1)).^2)/2).^(-lambda) ) );
end

norm(u(:,1) - u0_exact)
hold on; plot3(xgrid, zeros(size(xgrid)), u0_exact, '-or');


