clear all
close all    
clc
figure   


N = 30;
k = 1:N;
x = cos(((2*k -1)* pi)/(2*N));
y = rungef(x);
T= zeros(N,N);

for n = 0:N-1
    T(:,n+1) = cos(n .* acos(x));
end

c = T \ y';

xspace = linspace(-1,1,1000);
Tx = zeros(length(xspace), N);
for n = 0:N-1
    Tx(:,n+1) = cos(n * acos(xspace));
end

irungef = Tx * c;

a = dct(y, 'Type', 2);
p = ones(N,1);
p(1) = .5;
c2 = (a(:) .* p(:)) * 2 / N;

irungef_dct = Tx * c2;
norm(c - c2, Inf)





plot(xspace, rungef(xspace), 'k', 'LineWidth', 2); hold on;
plot(xspace, irungef_dct, 'r--', 'LineWidth', 1.5);
plot(xspace, irungef, 'b:', 'LineWidth', 1.5);
plot(x, y, 'go');

legend('True Runge', 'DCT Interpolant', 'Vandermonde Interpolant', 'Chebyshev Nodes');
xlabel('x'); ylabel('f(x)');
title('Chebyshev Interpolation of the Runge Function');
grid on;




function f = rungef(x)
f = 1 ./ (1+25.*x.^2);
end