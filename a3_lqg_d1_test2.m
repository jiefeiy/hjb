clear, clc
% implement LQG control problem for d=1
T = 1;
lambda = 1;
g = @(x) max(1-abs(x), 0);  % terminal condition 

% decide truncation boundary 
xmin = -6;
xmax = -xmin;
mmin = -1;
mmax = 1;

% discretization
n_time = 30;
dt = T/n_time;


knots1 = @(n) knots_CC(n,xmin,xmax,'nonprob');
knots2 = @(n) knots_CC(n,mmin,mmax,'nonprob');
S = create_sparse_grid(2,5,{knots1, knots2},@lev2knots_doubling); % grid
Sr=reduce_sparse_grid(S);

nq = 2^5;     % number of Hermite-Gauss knots 
wknots = @(n) knots_normal(n,0,1);
wS = create_sparse_grid(1, nq-1, wknots, @lev2knots_lin);

xq_grid = Sr.knots(1,:) + 2*sqrt(lambda)*Sr.knots(2,:)*dt + sqrt(2*dt)*wS.knots'; % nq-by-Sr.size matrix
m_star_on_xq = zeros(size(xq_grid));

%% dynamic programming: interpolate Q function
V_on_xq = g(xq_grid);
for n = n_time-1:-1:0
    Q_on_grid = Sr.knots(2,:).^2 + wS.weights*(V_on_xq);
    Q_fun = @(x,m) interpolate_on_sparse_grid(S,Sr,Q_on_grid,[x;m]);
    V_on_xq = zeros(size(xq_grid));
    for i = 1:nq
        for j = 1:Sr.size
            if xq_grid(i,j) > xmin && xq_grid(i,j) < xmax
                if n>0
                    [~,V_on_xq(i,j)] = fminbnd(@(m) Q_fun(xq_grid(i,j),m), mmin, mmax);
                else
                    [m_star_on_xq(i,j), V_on_xq(i,j)] = fminbnd(@(m) Q_fun(xq_grid(i,j),m), mmin, mmax);
                end
            else
                V_on_xq(i,j) = g(xq_grid(i,j));
            end
        end
    end
    n
end

%%
figure(1);
plot(xq_grid(:), V_on_xq(:), '.');

figure(2);
plot(xq_grid(:), m_star_on_xq(:), '.');

%% exact solution approximation 
% exact solution: u(0,x) = -1/lambda * log( E[ exp( - lambda*g(x + sqrt(2)*W_T) ) ] )
rng("default");
n_x = 50;
dx = (xmax - xmin)/n_x;
xgrid = (xmin:dx:xmax)';
M = 100000;
u0_exact = zeros(size(xgrid));
for j = 1:n_x+1
    u0_exact(j) = -1/lambda * log( mean( exp(-lambda*g(xgrid(j) + sqrt(2)*sqrt(T)*randn(M,1))) ) );
end

%%
figure(1);
hold on; plot(xgrid, u0_exact, '-r');

xlabel('x')
ylabel('$V_0(x)$', Interpreter='latex');
legend('naive interpolate Q','exact')

%% compute exact m
m_exact = zeros(n_x+1, 1);
for j = 2:n_x
    m_exact(j) = -sqrt(lambda)*(u0_exact(j+1) - u0_exact(j-1))/(2*dx);
end

figure(2);
hold on; plot(xgrid, m_exact, '-r');

xlabel('x')
ylabel('$m_0(x)$', Interpreter='latex');
legend('naive interpolate Q','exact')