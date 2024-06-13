clear, clc, close all
% use Feynman-Kac formula

% N = 512;
Nl = 32;
d = 2;
tau = 3;
% f = @(x) 1000*(x<=0.5) + 2000*(x>0.5);         % source of fluid 
f = @(x) ones(size(x));

%% KL expansion 
phi = @(l,x) sqrt(2)*cos(pi*l.*x);
lambda = @(l) (pi^2*l.^2 + tau.^2).^(-d);
sqrt_lambda_all = sqrt(lambda(1:Nl));

% theta_all = randn(1, Nl);                      % randomness comes from here
load('theta_all.mat');
a = @(x) exp( sum(theta_all .* sqrt_lambda_all .* phi((1:Nl), x), 2) );

phi_dx = @(l,x) -sqrt(2)*pi*l .* sin(pi*l.*x);
a_dx = @(x) sum(theta_all .* sqrt_lambda_all .* phi_dx((1:Nl), x), 2) .* a(x);

%%% plot one realization of permeability field a(x,theta)
% xgrid = linspace(0, 1, N+1)';
% plot(xgrid, a(xgrid), 'DisplayName', 'permeability');
% xlabel('x'); legend;
% figure(2);
% plot(xgrid, a_dx(xgrid), 'DisplayName', 'permeability derivative');
% xlabel('x'); legend;

%% generate paths of sde
xq = 2/3;     
M = 10000;
N_time = 20000;
T = 1;
dt = T / N_time;
dw_sample = randn(M, N_time) * sqrt(dt);
x_sample = zeros(M, N_time+1);
x_sample(:,1) = xq;
for k = 1:N_time
    x_sample(:, k+1) = x_sample(:,k) + a_dx(x_sample(:,k)) * dt ...
                        + sqrt(2*a(x_sample(:,k))) .* dw_sample(:,k);
end

% plot(0:N_time, x_sample(1:20, :))

k_first_exit = N_time * ones(M,1);
f_sample = zeros(size(x_sample));
for m = 1:M
    result = find(abs(x_sample(m,:) - 0.5)>0.5, 1);
    if ~isempty(result)
        k_first_exit(m) = result;
        f_sample(m, 1:result-2) = f(x_sample(m, 1:result-2));
    end
end
find(k_first_exit == N_time)

u_xq = sum(f_sample, "all")*dt/M

% save('theta_all.mat', "theta_all");

