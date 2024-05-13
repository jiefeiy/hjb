clear, clc
% numerical tests:
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
[X1, X2] = ndgrid(x1_grid, x2_grid);
u = sin(T + X1 + X2);

x0 = [5;6];
x0 = mod(x0, 2*pi);
sin(sum(x0))

%% matrix free implementation
for t = 1:n_time
    unew = zeros(nx+1, nx+1);
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
        
                laplace_u = (u(ip,j) + u(im,j)...
                     + u(i,jp) + u(i,jm) - 4*u(i,j));

                unew(i,j) = u(i,j) - dt/(2*2*h)*(u(ip,j) - u(im,j) + u(i,jp) - u(i,jm)) ...
                     + dt/(2*h*h)*(sigma_max^2*max(laplace_u,0) - sigma_min^2*max(-laplace_u,0)) ...
                     + dt*(sigma_min^2*max(u(i,j),0) - sigma_max^2*max(-u(i,j),0));
        end
    end
    u = unew;
end

uq = interpn(X1,X2,u,x0(1),x0(2));
uq

u_exact = sin(X1 + X2);
surf(X1, X2, u - u_exact);



    
