%% ASEN 3128
%Lab 2 - Simulate EOM Quadrotor
%Matthew Pabin
%Gavin O'Connell
%Mary Hanson
%Sruthi Bandla
%Code Created: 9/7/21

%% House Keeping
close all
clear
clc

%% Create Constants
g = 9.8;
m = 0.068;
r = 0.060;
k_m = 0.0024;
I_x = 0.000068;
I_y = 0.000092;
I_z = 0.000135;
nu = 0.001;
mu = 0.000002;

%% Initialize vectors
inertia = [I_x,0,0;0,I_y,0;0,0,I_z];
xi = [0;0;-10;0;0;0;0;0;0;0;0;0];

forces_motor = [(m*g)/4;(m*g)/4;(m*g)/4;(m*g)/4];

tspan = [0 10];
[t,state] = ode45( @(t,x) quadFun(t,x,g,m,inertia,forces_motor,k_m,nu,mu,r), tspan, xi);

%% Function Junction
function state_dot = quadFun(t,x,g,m,inertia,forces_motor,k_m,nu,mu,r)
iv = [inertia(1,1);inertia(2,2);inertia(3,3)];

moments_aero = -mu*sqrt(x(10)^2 + x(11)^2 + x(12)^2)*x(10:12);
forces_aero = -nu * sqrt(x(7)^2 + x(8)^2 + x(9)^2)* x(7:9);

moments_c = [-(r/sqrt(2)),-(r/sqrt(2)),(r/sqrt(2)),(r/sqrt(2));(r/sqrt(2)),-(r/sqrt(2)),-(r/sqrt(2)),(r/sqrt(2));...
    k_m,-k_m,k_m,-k_m] * forces_motor;
forces_c = [0;0;(-sum(forces_motor))];

position_dot = [(cos(x(5))*cos(x(6))) , (sin(x(4))*sin(x(5))*cos(x(6))) - (cos(x(4)) * sin(x(6))), (cos(x(4))*sin(x(5))* cos(x(6))) + (sin(x(4))*sin(x(6))); ...
    (cos(x(5))*sin(x(6))) , (sin(x(4))*sin(x(5))*sin(x(6))) + (cos(x(4))*cos(x(6))), (cos(x(4))*sin(x(5))*sin(x(6))) - (sin(x(4))*cos(x(6))); ...
    -sin(x(5)) , sin(x(4))*cos(x(5)) , cos(x(4))*cos(x(5))] * x(7:9);

angles_dot = [1 sin(x(4))*tan(x(5)) cos(x(4))*tan(x(5)); 0 cos(x(4)) -sin(x(4)); 0 sin(x(4))*sec(x(5)) cos(x(4))*sec(x(5))] * x(10:12);

vel_dot = [x(12)*x(8)-x(11)*x(9); x(10)*x(11)-x(12)*x(7); x(11)*x(7)-x(10)*x(8)] + g*[-sin(x(5)); cos(x(5))*sin(x(4)); cos(x(5))*cos(x(4))] + 1/m * forces_aero + 1/m * forces_c;

angvel_dot = [(iv(2) - iv(3))/iv(1) * x(11) * x(12);...
    (iv(3) - iv(1))/iv(2) * x(10) * x(12);...
    (iv(1) - iv(2))/iv(3) * x(10) * x(11);] + ...
    [1/(iv(1)) * moments_aero(1); 1/iv(2) * moments_aero(2); 1/iv(3) * moments_aero(3)] + ...
    [1/(iv(1)) * moments_c(1); 1/iv(2) * moments_c(2); 1/iv(3) * moments_c(3)];

state_dot = [position_dot;angles_dot;vel_dot;angvel_dot];
end
