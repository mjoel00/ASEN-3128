close all
clear
clc

%% Create Constants
global k_p
g = 9.8;
m = 0.068;
r = 0.060;
k_m = 0.0024;
k_p =-0.004;
I_x = 0.000068;
I_y = 0.000092;
I_z = 0.000135;
nu = 0.001;
mu = 0.000002;

%% Initialize vectors
inertia = [I_x,0,0;0,I_y,0;0,0,I_z];
tspan = [0 5];
x_initial = [0;0;-10; 0;0;0; 0;0;0; 0;0;0]; % [xe;ye;ze;phi;theta;psi;ue;ve;we;p;q;r]
A = [0,0,0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0,1,0,0,0;...
    0,0,0,0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,0,0,1;...
    0,0,0,0,-g,0,0,0,0,0,0,0;0,0,0,g,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;...
    0,0,0,-k2_lat/I_x,0,0,0,0,0,-k1_lat/I_x,0,0;0,0,0,0,-k2_long/Iy,0,0,0,0,0 ...
    ,-k2_long/Iy,0;0,0,0,0,0,0,0,0,0,0,0,0];
B = zeros(12);
B(9,9) = g;
C = eye(12);
D = zeros(12);
sys = ss(A,B,C,D);
[y,j,x] =initial(sys,x_initial);
forces_motor = ComputeMotorForces(m*g,0,0,0,r,k_m);

[t,state1] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor), tspan, x_initial);

j= t(1:21);
cont = [m*g,0,0,0]';

PlotAircraftSim(j,y,cont,'r')

function state_dot = quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor)
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


function motor_forces = ComputeMotorForces(Zc,Lc,Mc,Nc,R,km)

c = inv([-1,-1,-1,-1;-R/sqrt(2),-R/sqrt(2),R/sqrt(2),R/sqrt(2); ...
    R/sqrt(2),-R/sqrt(2),-R/sqrt(2),R/sqrt(2); km,-km,km,-km]);

motor_forces = c * [Zc; Lc; Mc; Nc];

end

function [] = PlotAircraftSim(time, aircraft_state_array, control_input_array, col)
%{
Description: PlotAirCraftSim plots information about the simulation whose
information is contained in aircraft_state_array and time.

Inputs: 
    time = time vector measured in [s].

    aircraft_state_array = a n x 12 state array. Each column represents one
    state of the aircraft over the simulation. In order:
        Columns 1:3 = X, Y, Z   [m]
        Columns 4:6 = φ, Θ, ψ   [rad]
        Columns 7:9 = u, v, w   [m/s]
        Columns 10:12 = p, q, r [rad/s]

    control_input_array = a n x 4 vector containing the control inputs to
    the four quadrotor motors over time. Each column represents the inputs
    to the corresponding motor.

    col = plotting format as a string. Dictates line and point style. 

Outputs: 
    Many, many figures. Note that these will not be saved to the workspace.
    
    
Methodology:
    Plot each figure independently, using hold on and hold off for clarity.
%} 

%% Position Figure
    figure(1)
% X
    hold on
    subplot(3,1,1)
    plot(time,(aircraft_state_array(:,1)),col);
    ylabel('x [m]');
    xlabel('Time [s]');
    hold off
% Y
    hold on
    subplot(3,1,2)      
    plot(time,(aircraft_state_array(:,2)),col);
    ylabel('y [m]');
    xlabel('Time [s]');
    hold off
% Z
    hold on
    subplot(3,1,3)      
    plot(time,(aircraft_state_array(:,3)),col);
    ylabel('z [m]');
    xlabel('Time [s]');
    hold off
    
%% Euler Angle Figure
    figure(2)
    hold on
% Roll
    hold on
    subplot(3,1,1)      
    plot(time,(aircraft_state_array(:,4)),col);
    ylabel('\phi [rad]');
    xlabel('Time [s]');
    hold off
% Pitch
    hold on
    subplot(3,1,2)      
    plot(time,(aircraft_state_array(:,5)),col);
    ylabel('\theta [rad]');
    xlabel('Time [s]');
    hold off
% Yaw
    hold on
    subplot(3,1,3)      
    plot(time,(aircraft_state_array(:,6)),col);
    ylabel('\psi [rad]');
    xlabel('Time [s]');
    hold off
    
%% Velocity Figure    
    figure(3)
% u 
    hold on
    subplot(3,1,1)      
    plot(time,(aircraft_state_array(:,7)),col);
    ylabel('u [m/s]');
    xlabel('Time [s]');
    hold off
% v 
    hold on
    subplot(3,1,2)      
    plot(time,(aircraft_state_array(:,8)),col);
    ylabel('v [m/s]');
    xlabel('Time [s]');
    hold off
% w 
    hold on
    subplot(3,1,3)
    plot(time,(aircraft_state_array(:,9)),col);
    ylabel('w [m/s]');
    xlabel('Time [s]');
    hold off
    
%% Angular Velocity Figure
    figure(4)
% p
    hold on
    subplot(3,1,1)
    plot(time,(aircraft_state_array(:,10)),col);
    ylabel('p [rad/s]');
    xlabel('Time [s]');
    hold off
% q
    hold on
    subplot(3,1,2)      
    plot(time,(aircraft_state_array(:,11)),col);
    ylabel('q [rad/s]');
    xlabel('Time [s]');
    hold off
% r
    hold on
    subplot(3,1,3)      
    plot(time,(aircraft_state_array(:,12)),col);
    ylabel('r [rad/s]');
    xlabel('Time [s]');
    hold off
    
%% Control Input Figure
    figure(5)
% First motor
    hold on
    subplot(4,1,1)
    plot(time,control_input_array(:,1),col);
    ylabel('F_1 [N]');
    xlabel('Time [s]');
    hold off
% Second motor
    hold on
    subplot(4,1,2)
    plot(time,control_input_array(:,2),col);
    ylabel('F_2 [N]');
    xlabel('Time [s]');
    hold off
% Third motor
    hold on
    subplot(4,1,3)
    plot(time,control_input_array(:,3),col);
    ylabel('F_3 [N]');
    xlabel('Time [s]');
    hold off
% Fourth motor
    hold on
    subplot(4,1,4)
    plot(time,control_input_array(:,4),col);
    ylabel('F_4 [N]');
    xlabel('Time [s]');
    hold off
    
%% 3D Flight Path
    figure(6)
    %Extract X, Y and Z position data from aircraft_state_array
    xData = aircraft_state_array(:,1);
    yData = aircraft_state_array(:,2);
    zData = -1*aircraft_state_array(:,3); %Correct from body frame 
    %Plotting
    hold on
    grid on
    view(-13,12)
    % Do fancy colors ;)
    scatter3(xData,yData,zData,3,-zData)
    colormap(jet);
    %This next part gets tricky, I want to fill the above with the
    %corresponding color of the trace (ie, the color information contained
    %in col). So...
    %...create a dummy plot to extract color info from col
    A = plot(0,0,col);
    colColor = get(A,'Color');
    % Add points to the starting and ending positions
    B = plot3(xData(end),yData(end),zData(end),'rp','MarkerFaceColor',colColor,'linewidth',1.5,'markersize',15);
    C = plot3(xData(1),yData(1),zData(1),'gp','MarkerFaceColor',colColor,'linewidth',1.5,'markersize',15);
    %Add the fill colors
    B.MarkerFaceColor = colColor;
    C.MarkerFaceColor = colColor;
    % Add labels
    ylabel('Y Position [m]');
    xlabel('X Position [m]');
    zlabel('Elevation [m]');
    hold off
end
