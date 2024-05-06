clear, clc
% numerical tests:
% [1] Forsyth, Labahn. (2007) Numerical methods for controlled
% Hamilton-Jacobi-Bellman PDEs in finance. 
%
% Example 2.3: American option using Policy iteration in Section 6.2
% \min( -V_t - {\sigma^2 S^2/2 V_{SS} + rS V_S - rV}, V - V^* ) = 0.
% 
% This is equivalent to a penalized control problem 
% V_t + 
% \sup_{\mu\in\{0, 1\}} {\sigma^2 S^2/2 V_{SS} + rS V_S - rV + \mu(V^* - V)/\varepsilon} = 0.
%
% Let L^Q V := a(S, \tau, Q) V_SS + b(S, \tau, Q)V_s - c(S, \tau, Q)V.
% The HJB has the form of 
% V_\tau = sup_Q {L^Q V + d(S, \tau, Q)}. 
% 
% - a(S, \tau, Q) = \sigma^2 S^2/2,
% - b(S, \tau, Q) = rS, 
% - c(S, \tau, Q) = r + \mu/\varepsilon
% - d(S, \tau, Q) = \mu/varepsilon
%
% Aim to solve the problem on [0,T]\times[0, S_{max}]

%% American put parameters
p.price = 36;
p.volatility = 0.2;
p.rate = 0.06;
p.expiration = 1;
p.strike = 40;
eps = 0.01;   % penalty level

%% determine the boundary condition
price_max = 100;
V_star = @(s) max(p.strike - s, 0);

%% discretization parameters
n_time = 1000;
n_price = 500;
dt = p.expiration / n_time;
h = price_max / n_price;
price_grid = (0:h:price_max)';

%% grid and matrix
sigma = p.volatility;
r = p.rate;
alpha = sigma^2/(2*h^2) * price_grid.^2 - r/(2*h) * price_grid;
beta = sigma^2/(2*h^2) * price_grid.^2 + r/(2*h) * price_grid;

%% solve the algebraic discrete equations
tol = 1e-3;
max_iter = 5; num_iter_all = zeros(n_time, 1);

V = zeros(n_price+1, n_time+1);
V(:, 1) = V_star(price_grid);  % initial value
for n = 1:n_time
    Vtemp = V(:, n);
    k=1; criteria = 1;
    while k <= max_iter && criteria > tol     % policy iteration
        %% solve maximization problem
        %%% Case 1: mu = 0
        mu = zeros(n_price+1, 1);
        AQ = diag(-alpha - beta - r - mu/eps) + diag(alpha(2:end), -1) + diag(beta(1:end-1), 1);
        AQ(end, :) = 0;               % last row = 0 for boundary condition
        DQ = mu/eps .* V_star(price_grid);
        obj_0 = (AQ * Vtemp + DQ);
        %%% Case 2: mu = 1
        mu = ones(n_price+1, 1);
        AQ = diag(-alpha - beta - r - mu/eps) + diag(alpha(2:end), -1) + diag(beta(1:end-1), 1);
        AQ(end, :) = 0;               % last row = 0 for boundary condition
        DQ = mu/eps .* V_star(price_grid);
        obj_1 = (AQ * Vtemp + DQ);
    
        idx = (obj_1 > obj_0);                       % compare Case 1 and 2
        mu = zeros(n_price+1, 1); mu(idx) = 1;       % find optimal control
        AQ = diag(-alpha - beta - r - mu/eps) + diag(alpha(2:end), -1) + diag(beta(1:end-1), 1);
        AQ(end, :) = 0;               % last row = 0 for boundary condition
        DQ = mu/eps .* V_star(price_grid);
        %% compute stopping criteria
        Vprep = Vtemp;
        Vtemp = (eye(n_price+1) - dt*AQ) \ (V(:, n) + dt*DQ);
        criteria = max(abs(Vtemp - Vprep));
        k=k+1;
    end
    num_iter_all(n) = k;                  % record the number of iterations
    V(:, n+1) = Vtemp;
end

%% display the results
plot(0:h:price_max, V(:,end), '.-k');
hold on; plot(0:h:price_max, V(:, 1), '.-r');

% [X,Y] = meshgrid(0:dt:p.expiration, 0:h:price_max);
% surf(X, Y, V);

amer = interp1(price_grid, V(:,end), p.price)

[~,euro_ref] = blsprice(p.price, p.strike, p.rate, p.expiration, p.volatility)






