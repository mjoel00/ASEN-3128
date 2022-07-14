% Matthew Pabin
% ASEN 3128
% Simulate Equations of Motion 2
% 8/24/21
 
clear
clc
 
%% Variables
 
rho = 0.961; %kg/m^3 ambient air density
Cd = 0.6; %coefficient of drag
A = 2*pi*0.015^2; %half sphere of golf ball, front surface area
g = 9.81; %9.81m/s^2 acceleration downwards (gravity)
golf0 = [0,0,0,0,20,20]'; %coordinate system x-north, y-east, z-up for graphing purposes
tspan = [0 5]; %ODE45 integration time
x = cell(1,100);
t = cell(1,100);
 
%% Generating Cells of random and set variables for different ODE45 loops
 
wind = cell(0,100);
wind_2 = cell(0,100);
for i = 1:100
    wind{i} = [0 0 0]; %1x100 cells containing 0 windspeed
    wind_2{i} = [0.05*randi([-100,100]),0.05*randi([-100,100]),0]; %1x100 cell with independent random x and y velocity between -5m/s and 5m/s
end
 
m = cell(0,100);
m_2 = cell(0,100);
for i = 1:100
    m{i} = 0.6; %1x100 cells containing 0.6g mass
    m_2{i} = 0.1+0.01*randi([0,490]); %1x100 cells with random masses between 0.1 and 5g
end
 
%% Stopping function call for when projectile hits x = 0 for each ODE45 set
opt = odeset('Events',@stop);
 
%% ODE call for wind variation
for i = 1:100
    [t{i},x{i}] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m{i},g,wind_2{i}),tspan,golf0,opt);
end
 
%3D plot of flight path affected by wind variations between 5m/s in each direction
figure 
plot3(x{1}(:,1),x{1}(:,2),x{1}(:,3));
title('3D Flight Path with Random Wind Variation (-5m/s to 5m/s on x and y)');
hold on
grid on
xlabel('x (meters)')
ylabel(' y (meters)')
zlabel('z (meters)');
for i=1:99
   plot3(x{i}(:,1),x{i}(:,2),x{i}(:,3)) 
end
hold off
 
%% ODE call for mass variation
for i = 1:100
    [t{i},x{i}] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m_2{i},g,wind{i}),tspan,golf0,opt);
end
 
%2D plot of flight path affected by mass variations of the ball
figure 
plot(x{1}(:,2),x{1}(:,3));
title('2D Flight path with Random Mass Variation (0.1g to 5g)');
hold on
grid on
xlabel('Distance')
ylabel('Height')
 
for i=1:99
   plot(x{i}(:,2),x{i}(:,3)) 
end
hold off
%% Landing Points
landingpts = zeros(100,3);
for i = 1:100
    landingpts(i,:) = x{1,i}(end,1:3);
end 

 
%% ODE Function
function x = objectEOM(t,x,rho,Cd,A,m,g,wind)
    
    %applying windspeed to x and y directions
    inertial_airspeed(1) = x(4) + wind(1);
    inertial_airspeed(2) = x(5) + wind(2);
    %calculating drag in each velocity component
    F(1) = -A * Cd * 0.5 * rho * inertial_airspeed(1)^2; 
    F(2) = -A * Cd * 0.5 * rho * inertial_airspeed(2)^2;
    F(3) = -A * Cd * 0.5 * rho * x(6)^2;
    a = F / m;          %turn force to acceleration
    a(3) = a(3) - g;    %apply gravity acceleration to the z direction
    
    x(1) = inertial_airspeed(1);
    x(2) = inertial_airspeed(2);
    x(3) = x(6);
    x(4) = a(1);
    x(5) = a(2);
    x(6) = a(3);
end
 
%% Stopping function
 
function [value, isterminal, direction] = stop(t, x)
    value = x(3);
    isterminal = 1;
    direction = -1;
end
