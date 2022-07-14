%ASEN 3128 Lab IV

%Matlab Maintenance
clear, clc, close all

%Constants

c.Ix = 5.8*10^-5;     %[kg*m^2] Inertia in the body x
c.Iy = 7.2*10^-5;     %[kg*m^2] Inertia in the body y
c.m = 0.068;          %[kg] Mass
c.g = 9.81;           %[m/s^2]
sf = 0.1;           % Scale Factor; relates tau1, tau2
Tau_1 = 0.5;        %[s] Time constant (given)
Tau_2 = Tau_1*sf;   %[s] Second Time constant (picked by yours truely)



%Math
lambda_1 = -1/Tau_1;    % Eigenvalue for tau1
lambda_2 = -1/Tau_2;    % Eigenvalue for tau2

k2_lat = c.Ix*(lambda_2*lambda_1);
k1_lat = c.Ix*(-lambda_2-lambda_1);

k2_long = c.Iy*(lambda_2*lambda_1);
k1_long = c.Iy*(-lambda_2-lambda_1);


%% Linearized stuff
fig_a = [1,2,3,4,5,6]; % figure a number vector
fig_b = [1,2,3,4,5,6] + 6; % figure b number vector
fig_c = [1,2,3,4,5,6] + 2*6; % figure c number vector
fig_d = [1,2,3,4,5,6] + 3*6; % figure d number vector
fig_e = [1,2,3,4,5,6] + 4*6; % figure e number vector
fig_f = [1,2,3,4,5,6] + 5*6; % figure f number vector
control = '0.1 rad/s pitch rate';
col='c-';
col2='m.';


A = [0,0,0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0,1,0,0,0;...
    0,0,0,0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,0,0,1;...
    0,0,0,0,-c.g,0,0,0,0,0,0,0;0,0,0,c.g,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;...
    0,0,0,-k2_lat/c.Ix,0,0,0,0,0,-k1_lat/c.Ix,0,0;0,0,0,0,-k2_long/c.Iy,0,0,0,0,0 ...
    ,-k2_long/c.Iy,0;0,0,0,0,0,0,0,0,0,0,0,0];

B = zeros(12);
B(9,9) = c.g;
C = eye(12);
D = zeros(12);
x0=[0;0;0;deg2rad(5);0;0;0;0;0;0;0;0];
system=ss(A,B,C,D);
[x,t,y]=initial(system,x0,5);
PlotAircraftSim(t,y,control,col,fig_a)

    

    
%% Non-linear portion
g = 9.81; %[m/s^2]
m = 0.068; %[kg]
R = 0.06; %[m]
km = 0.0024; %[N*m/N]
Ix = 5.8*10^(-5); %[kg*m^2]
Iy = 7.2*10^(-5); %[kg*m^2]
Iz = 1*10^(-4); %[kg*m^2]
vCoef = 1*10^(-3); %[N/(m/s)^2]
mu = 2*10^(-6); %[N*m/(rad/s)^2]

C = {g m R km Ix Iy Iz vCoef mu k1_lat k2_lat k1_long k2_long};
[t, s] = ode45(@(t,x) odeSimControlled(t,x,C),[0 5],x0);


PlotAircraftSim(t,s,control,col2,fig_a)