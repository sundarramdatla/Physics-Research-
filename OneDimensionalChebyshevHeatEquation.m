clear all
clc

alpha = 1;
S = 1;
N = 32;
tspan = [0,2];

k = 0:N-1;
x = cos(pi * k / (N-1))';

u0 = zeros(N,1);
u0(1) = 1;
u0(end) = 0.0;

[t,U] = ode45(@(t,u) my_pde(u, alpha, S, N), tspan, u0);

num_snapshots = 5;
t_indices = round(linspace(1, length(t), num_snapshots));

figure;
hold on;
colors = lines(num_snapshots); 

for i = 1:num_snapshots
    idx = t_indices(i);
    plot(x, U(idx, :), 'LineWidth', 2, 'Color', colors(i,:), ...
        'DisplayName', sprintf('t = %.2f', t(idx)));
end

xlabel('x'); ylabel('u(x,t)');
title('Heat Equation Snapshots at Different Times');
legend show;
grid on;

function dudt = my_pde(u, alpha, S, N)
p = ones(N,1);
p(1) = 0.5;
p(end) = .5;

a = dct(u, 'Type', 1);
c = a .* p / (N-1);
c_dd = zeros(N,1);
    for k = 0:N-1
        c_dd(k+1) = -k^2 * c(k+1);  
    end

a_dd = (N-1) * c_dd ./ p;

u_xx = idct(a_dd, 'Type', 1);

dudt = alpha * u_xx + S;
dudt(1) = 0;
dudt(end) = 0;




end





disp(U(end, :))