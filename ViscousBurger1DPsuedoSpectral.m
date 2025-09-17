clear all
clc

t0 = 0;
tf = 50;
t = linspace(t0,tf,100);
v = .1;
x0 = 0;
xf = 10;
nx = 32;
x = linspace(x0,xf,nx + 1);
x(end) = [];
dx = x(2) - x(1);
L = xf - x0;

u0 = cos(2*pi*x / L) + 2;
u0 = u0(:);


options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', 0.1);
[t,u] = ode45(@my_pde, t, u0, options, L, v);


figure;
for i = 1:length(t)
    plot(x, u(i,:), 'LineWidth', 2)
    xlabel('x'); ylabel('u(x,t)');
    title(['Time t = ', num2str(t(i))])
    axis([x0 xf 0 4])
    grid on
    drawnow
end

function ut = my_pde(t,u, L, v)
n = length(u);
k = (2 * pi / L) * [0:(n/2 - 1), -n/2:-1]';
u_x = real(ifft(1i * k .* fft(u)));
u_xx = real(ifft(-k.^2 .* fft(u)));
ut = - u .* u_x + v * u_xx;
end



