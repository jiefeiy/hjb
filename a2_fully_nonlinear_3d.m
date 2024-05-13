clear, clc
% numerical tests:
% [1] W Guo, J Zhang, J Zhuo. (2015) A monotone scheme for high-dimensional
% fully nonlinear pdes.
%
% test with example 1 in Section 6 Numerical examples. 
% finite difference method with naive implementation.

%% parameter setting 
sigma_min = 1;
sigma_max = sqrt(2);
T = 0.5;

%% discritization and terminal condition
n_time = 200;
nx = 50;
dt = T / n_time;
h = 2*pi / nx;
x1_grid = linspace(0, 2*pi, nx+1);
x2_grid = linspace(0, 2*pi, nx+1);
x3_grid = linspace(0, 2*pi, nx+1);
[X1, X2, X3] = ndgrid(x1_grid, x2_grid, x3_grid);
u = sin(T + X1 + X2 + X3);

x0 = [5;6;7];
% id_grid = [20 20 5];
% x0 = [x1_grid(id_grid(1)); x2_grid(id_grid(2)); x3_grid(id_grid(3))];
x0 = mod(x0, 2*pi);
sin(sum(x0))

%% matrix free implementation
for t = 1:n_time
    unew = zeros(nx+1, nx+1, nx+1);
    for i = 1:nx+1
        if i == 1
            im = nx; ip = i+1;
        elseif i == nx+1
            im = i-1; ip = 2;
        else
            im = i-1; ip = i+1;
        end
        for j = 1:nx+1
                if j == 1
                    jm = nx; jp = j+1;
                elseif j == nx+1
                    jm = j-1; jp = 2;
                else
                    jm = j-1; jp = j+1;
                end
                for k = 1:nx+1
                        if k == 1
                            km = nx; kp = k+1;
                        elseif k == nx+1
                            km = k-1; kp = 2;
                        else
                            km = k-1; kp = k+1;
                        end
                        
                        laplace_u = (u(ip,j,k) + u(im,j,k)...
                             + u(i,jp,k) + u(i,jm,k) + u(i,j,kp) + u(i,j,km) - 6*u(i,j,k));
 
                        unew(i,j,k) = u(i,j,k) - dt/(3*2*h)*(u(ip,j,k) - u(im,j,k) + u(i,jp,k) - u(i,jm,k)...
                             + u(i,j,kp) - u(i,j,km)) + dt/(2*h*h)*(sigma_max^2*max(laplace_u,0) - sigma_min^2*max(-laplace_u,0)) ...
                             + dt*3/2*(sigma_min^2*max(u(i,j,k),0) - sigma_max^2*max(-u(i,j,k),0));
                end
        end
    end
    u = unew;
end

% u(id_grid(1), id_grid(2), id_grid(3))

uq = interpn(X1,X2,X3,u,x0(1),x0(2),x0(3));
uq



    
