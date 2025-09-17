
g = 9.81;
m = .27; 
R = 0.02; 
N = m*[0; 0; g];
I = [2/5 * m * R^2 0 0 ; 0 2/5 * m * R^2 0; 0 0 2/5 * m * R^2]; %inertia tensor
u_k = 0.4; 
e = 0.8; 
p = 1.2;
A = .001256; 
C = 0.47; 



F = [50;0;0]; 
F_x = F(1);
F_y = F(2);
F_z = F(3);

w0 = [0; 0; -10]; 
vx0 = (F_x * .1) / m;  
vy0 = (F_y * .1) / m; 
vz0 = (F_z *.1)/ m;
s0 = [0, vx0, 0, vy0, 5, vz0]; 
t_span = linspace(0,10,500); 



function dsdt = ode_func(s,m, g, p, A, C, w)

x1 = s(1);
x2 = s(2);
y1 = s(3);
y2 = s(4);
z1 = s(5);
z2 = s(6);
v = [x2;y2;z2];
v_mag = norm(v);

if v_mag ~= 0 
    F_d = .5*A*C*p*(v_mag^2) * (-v/v_mag);
else
    F_d = [0;0;0]
end
v_mag;
Re = 2*.02 * 1.2 * v_mag / (1.8 * 10^-5)
F_d_x = F_d(1);
F_d_y = F_d(2);
F_d_z = F_d(3);

F_m = C*cross(w, F_d); 
F_m_x = F_m(1);
F_m_y = F_m(2);
F_m_z = F_m(3);

dx1dt = x2;
dx2dt = F_m_x + F_d_x;
dy1dt = y2;
dy2dt = F_m_y + F_d_y;
dz1dt = z2;
dz2dt =  F_m_z + F_d_z -g;
dsdt = [dx1dt;dx2dt; dy1dt; dy2dt;dz1dt;dz2dt];
end



function [check, isterminal, direction] = groundFunc(t,s)
check = s(5);
isterminal = 1;
direction = -1;

end

options = odeset('Events', @groundFunc, 'Refine', 10);



[t, results_s] = ode45(@(t,s) ode_func(s, m, g, p, A, C,w0), t_span, s0, options);
tplot = t;
splot = results_s;
wplot = w0';



function dwdt = ode_func2(w, T, I) 
I_inv = inv(I);
dwdt = I_inv * T;
end

function ds2dt = ode_func3(s, u_k, g) 
x1 = s(1);
x2 = s(2);
y1 = s(3);
y2 = s(4);
z1 = s(5);
z2 = s(6);
dx1dt = x2;
dx2dt = - (u_k * g * x2)/ (sqrt(x2^2 + y2^2));
dy1dt = y2;
dy2dt = - (u_k * g * y2)/ (sqrt(x2^2 + y2^2));
dz1dt = z2;
dz2dt = -g;
ds2dt = [dx1dt; dx2dt; dy1dt; dy2dt; dz1dt; dz2dt];
end



for i = 1:10
t_start = tplot(end);
t_end = t_start + 10;
s_last = splot(end, :);
F_k = -u_k * m * g * ([s_last(2);s_last(4); 0]/ norm([s_last(2); s_last(4); 0])); 
 
[t, results_s_f] = ode45(@(t,s_f) ode_func3(s_f, u_k, g), [t_start t_start + .01], [s_last(1), s_last(2), s_last(3), s_last(4), 0, 0]);
splot = [splot; results_s_f]; 
s_next = [results_s_f(end, 1), results_s_f(end, 2), results_s_f(end, 3), results_s_f(end, 4), 0, -e * s_last(6)];




r = -R * [0;0;1]; 
T = cross(F_k, r); 
[t, results_w] = ode45(@(t,w) ode_func2(w, T, I), [t_start t_start + 0.01], wplot(end, :)');
wplot = [wplot; results_w(end, :)];


[t, results_s] = ode45(@(t,s) ode_func(s, m, g, p, A, C, wplot(end, :)'), [t_start t_end], s_next, options);
tplot = [tplot; t];
splot = [splot;results_s];


end

figure
comet3(splot(:,1), splot(:,3), splot(:,5))

figure
comet(splot(:,1), splot(:,5))
title('Ping pong ball in 3D Space')


ylabel('y')
xlabel('x')
zlabel('z')
grid on

figure
title('Top Down View of Trajectory')
plot(splot(:,1), splot(:,3))

