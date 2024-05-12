clear, clc
% numerical tests:
% [1] W Guo, J Zhang, J Zhuo. (2015) A monotone scheme for high-dimensional
% fully nonlinear pdes.
%
% test with example 1 in Section 6 Numerical examples. 

%% parameter setting 
sigma_min = 1;
sigma_max = sqrt(2);
T = 0.5;
x0 = [5 6 7]';
x0 = mod(x0, 2*pi);
sin(sum(x0))

%% discritization and terminal condition
n_time = 2000;
nx = 30;
dt = T / n_time;
h = 2*pi / nx;
x1_grid = linspace(0, 2*pi, nx+1);
x2_grid = linspace(0, 2*pi, nx+1);
x3_grid = linspace(0, 2*pi, nx+1);
[X1, X2, X3] = ndgrid(x1_grid, x2_grid, x3_grid);
u = sin(T + X1 + X2 + X3);

%% matrix free implementation
sigma = sigma_max;
for t = 1:n_time
    unew = zeros(nx+1, nx+1, nx+1);
    for i = 1:nx+1
        if i == 1
            im = nx+1; ip = i+1;
        elseif i == nx+1
            im = i-1; ip = 1;
        else
            im = i-1; ip = i+1;
        end
        for j = 1:nx+1
                if j == 1
                    jm = nx+1; jp = j+1;
                elseif j == nx+1
                    jm = j-1; jp = 1;
                else
                    jm = j-1; jp = j+1;
                end
                for k = 1:nx+1
                        if k == 1
                            km = nx+1; kp = k+1;
                        elseif k == nx+1
                            km = k-1; kp = 1;
                        else
                            km = k-1; kp = k+1;
                        end
    
                         unew(i,j,k) = u(i,j,k) + dt/(3*2*h)*(u(ip,j,k) - u(im,j,k) + u(i,jp,k) - u(i,jm,k)...
                             + u(i,j,kp) - u(i,j,km)) + dt/2*sigma^2/(h^2)*(u(ip,j,k) + u(im,j,k)...
                             + u(i,jp,k) + u(i,jm,k) + u(i,j,kp) + u(i,j,km) - 6*u(i,j,k)) ...
                             - dt*3/2*sigma^2 * u(i,j,k);
                end
        end
    end
    u = unew;
end

X1t = pagetranspose(X1);
X2t = pagetranspose(X2);
X3t = pagetranspose(X3); % convert to meshgrid
ut = pagetranspose(u);
Vq = interp3(X1t,X2t,X3t,ut,x0(1),x0(2),x0(3));
Vq


% [X1m,X2m,X3m] = meshgrid(x1_grid, x2_grid, x3_grid);
% isequal(X1t,X1m) & isequal(X2t,X2m) & isequal(X3t,X3m)

% %% inner grid
% i = 2:nx; j=i; k=i;
% sigma = sigma_min;
% unew = u;
% unew(i,j,k) = u(i,j,k) + dt/(3*2*h)*(u(i+1,j,k) - u(i-1,j,k) + u(i,j+1,k) - u(i,j-1,k) + u(i,j,k+1) - u(i,j,k-1)) ...
%     + dt/2*sigma^2/(h^2)*(u(i+1,j,k) + u(i-1,j,k) + u(i,j+1,k) + u(i,j-1,k) + u(i,j,k+1) + u(i,j,k-1) - 6*u(i,j,k)) ...
%     - dt*3/2*sigma^2 * u(i,j,k);
% 
% %% periodic boundary condition
% unew(1,j,k) = u(1,j,k) + dt/(3*2*h)*(u(2,j,k) - u(nx+1,j,k) + u(1,j+1,k) - u(1,j-1,k) + u(1,j,k+1) - u(1,j,k-1)) ...
%     + dt/2*sigma^2/(h^2)*(u(2,j,k) + u(nx+1,j,k) + u(1,j+1,k) + u(1,j-1,k) + u(1,j,k+1) + u(1,j,k-1) - 6*u(1,j,k)) ...
%     - dt*3/2*sigma^2 * u(1,j,k);
% unew(nx+1,j,k) = u(nx+1,j,k) + dt/(3*2*h)*(u(1,j,k) - u(nx,j,k) + u(nx+1,j+1,k) - u(nx+1,j-1,k) + u(nx+1,j,k+1) - u(nx+1,j,k-1)) ...
%     + dt/2*sigma^2/(h^2)*(u(1,j,k) + u(nx,j,k) + u(nx+1,j+1,k) + u(nx+1,j-1,k) + u(nx+1,j,k+1) + u(nx+1,j,k-1) - 6*u(nx+1,j,k)) ...
%     - dt*3/2*sigma^2 * u(nx+1,j,k);
% 
% unew(i,1,k) = u(i,1,k) + dt/(3*2*h)*(u(i+1,1,k) - u(i-1,1,k) + u(i,2,k) - u(i,nx+1,k) + u(i,1,k+1) - u(i,1,k-1)) ...
%     + dt/2*sigma^2/(h^2)*(u(i+1,1,k) + u(i-1,1,k) + u(i,2,k) + u(i,nx+1,k) + u(i,1,k+1) + u(i,1,k-1) - 6*u(i,1,k)) ...
%     - dt*3/2*sigma^2 * u(i,1,k);
% unew(i,nx+1,k) = u(i,nx+1,k) + dt/(3*2*h)*(u(i+1,nx+1,k) - u(i-1,nx+1,k) + u(i,1,k) - u(i,nx,k) + u(i,nx+1,k+1) - u(i,nx+1,k-1)) ...
%     + dt/2*sigma^2/(h^2)*(u(i+1,nx+1,k) + u(i-1,nx+1,k) + u(i,1,k) + u(i,nx,k) + u(i,nx+1,k+1) + u(i,nx+1,k-1) - 6*u(i,nx+1,k)) ...
%     - dt*3/2*sigma^2 * u(i,nx+1,k);
% 
% unew(i,j,1) = u(i,j,1) + dt/(3*2*h)*(u(i+1,j,1) - u(i-1,j,1) + u(i,j+1,1) - u(i,j-1,1) + u(i,j,2) - u(i,j,nx+1)) ...
%     + dt/2*sigma^2/(h^2)*(u(i+1,j,1) + u(i-1,j,1) + u(i,j+1,1) + u(i,j-1,1) + u(i,j,2) + u(i,j,nx+1) - 6*u(i,j,1)) ...
%     - dt*3/2*sigma^2 * u(i,j,1);
% unew(i,j,nx+1) = u(i,j,nx+1) + dt/(3*2*h)*(u(i+1,j,nx+1) - u(i-1,j,nx+1) + u(i,j+1,nx+1) - u(i,j-1,nx+1) + u(i,j,1) - u(i,j,nx)) ...
%     + dt/2*sigma^2/(h^2)*(u(i+1,j,nx+1) + u(i-1,j,nx+1) + u(i,j+1,nx+1) + u(i,j-1,nx+1) + u(i,j,1) + u(i,j,nx) - 6*u(i,j,nx+1)) ...
%     - dt*3/2*sigma^2 * u(i,j,nx+1);



    
