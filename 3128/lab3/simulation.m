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
    0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0];
B = zeros(12);
B(9,9) = g;
C = eye(12);
D = zeros(12);
sys = ss(A,B,C,D);
[y,j,x] =initial(sys,x_initial);
%j = j(1:2001);
%y = y(1:2001); 
L_c = x_initial(10) *  k_p;
M_c = x_initial(11) * k_p;
N_c = x_initial(12)* k_p;
forces_motor = ComputeMotorForces(m*g,L_c,M_c,N_c,r,k_m);
%forces_motor = [-m*g/4;-m*g/4;-m*g/4;-m*g/4];
[t,state1] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor), tspan, x_initial);
F_c =  [x_initial(10),x_initial(11),x_initial(12)].*k_p;
PlotAircraftSim(t,state1, "r")
hold on
PlotAircraftSim(j,y, "b")


%% Function Junction

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


function state_dot = quadFunQ2(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor)
global k_p
iv = [inertia(1,1);inertia(2,2);inertia(3,3)];

%L_c = x(10) * k_p;
%M_c = x(11) * k_p;
%N_c = x(12) * k_p;
moments_aero = -mu*sqrt(x(10)^2 + x(11)^2 + x(12)^2)*x(10:12);

forces_aero = -nu * sqrt(x(7)^2 + x(8)^2 + x(9)^2)* x(7:9);

moments_c = [-(r/sqrt(2)),-(r/sqrt(2)),(r/sqrt(2)),(r/sqrt(2));(r/sqrt(2)),-(r/sqrt(2)),-(r/sqrt(2)),(r/sqrt(2));...]
    k_m,-k_m,k_m,-k_m] * forces_motor;
forces_c = [0;0;(-sum(forces_motor))];

position_dot = [(cos(x(5))*cos(x(6))) , (sin(x(4))*sin(x(5))*cos(x(6))) - (cos(x(4)) * sin(x(6))), (cos(x(4))*sin(x(5))* cos(x(6))) + (sin(x(4))*sin(x(6))); ...
    (cos(x(5))*sin(x(6))) , (sin(x(4))*sin(x(5))*sin(x(6))) + (cos(x(4))*cos(x(6))), (cos(x(4))*sin(x(5))*sin(x(6))) - (sin(x(4))*cos(x(6))); ...
    -sin(x(5)) , sin(x(4))*cos(x(5)) , cos(x(4))*cos(x(5))] * x(7:9);

angles_dot = [1 sin(x(4))*tan(x(5)) cos(x(4))*tan(x(5)); 0 cos(x(4)) -sin(x(4)); 0 sin(x(4))*sec(x(5)) cos(x(4))*sec(x(5))] * x(10:12);

vel_dot = [x(12)*x(8)-x(11)*x(9); x(10)*x(11)-x(12)*x(7); x(11)*x(7)-x(10)*x(8)] + g*[-sin(x(5)); cos(x(5))*sin(x(4)); cos(x(5))*cos(x(4))] + 1/m * forces_aero ;

angvel_dot = [(iv(2) - iv(3))/iv(1) * x(11) * x(12);...
    (iv(3) - iv(1))/iv(2) * x(10) * x(12);...
    (iv(1) - iv(2))/iv(3) * x(10) * x(11);] + ...
    [1/(iv(1)) * moments_aero(1); 1/iv(2) * moments_aero(2); 1/iv(3) * moments_aero(3)] + ...
    [1/(iv(1)) * moments_c(1); 1/(iv(2)) * moments_c(2); 1/(iv(3)) * moments_c(3)] ;

state_dot = [position_dot;angles_dot;vel_dot;angvel_dot];
end


function [x,t] = PlotAircraftSim(time,aircraft_state_array,col)
% enter 'k' for col, will incorporate at a later time

% subplots
%creating arrays of labels for use
titles = ["X Postion","Y Postion","Z Postion","Euler Angle Phi","Euler Angles Theta","Euler Angle Psi",...
    "X Velocity","Y Velocity","Z Velocity","Angular Velocity P","Angular Velocity Q","Angular Velocity R"];
xlabels = "Time";
ylabels = [("[m]"),("[m]"),("[m]"),("[rad]"),("[rad]"),("[rad]"),("[m/s]"),...
    ("[m/s]"),("[m/s]"),("[rad/s]"),("[rad/s]"),("[rad/s]")];


for a = 1:3:length(aircraft_state_array(1,:))
    figure(a);
    subplot(3,1,1);
    plot(time,aircraft_state_array(:,a),'-.');
    xlim([0 5])
    hold on;
    xlabel(xlabels);
    ylabel(ylabels(a));
    title(titles(a));
    
    subplot(3,1,2);
    plot(time,aircraft_state_array(:,a+1),'-.');
    xlim([0 5])
    hold on;
    xlabel(xlabels);
    ylabel(ylabels{a+1});
    
    
    subplot(3,1,3);
    plot(time,aircraft_state_array(:,a+2),'-.');
    xlim([0 5])
    hold on;
    xlabel(xlabels);
    ylabel(ylabels{a+2});
    
end
%% input array left commented out until we can generate the array
% figure(4)
% xlim([0 5])
% subplot(4,1,3);
% plot(time,control_input_array(:,1));
% hold on;
% xlabel(xlabels);
% ylabel("Input 1");
% title("Input 1");
% 
% subplot(4,1,2);
% xlim([0 5])
% plot(time,control_input_array(:,2));
% hold on;
% xlabel(xlabels);
% ylabel("Input 2");
% title("Input 2");
% 
% subplot(3,1,3);
% xlim([0 5])
% plot(time,control_input_array(:,3));
% hold on;
% xlabel(xlabels);
% ylabel("Input 3");
% title("Input 3");
% 
% subplot(3,1,4);
% xlim([0 5])
% plot(time,control_input_array(:,4));
% hold on;
% xlabel(xlabels);
% ylabel("Input 4");
% title("Input 4");

figure(5)
xlim([0 5])
plot3(aircraft_state_array(:,1),aircraft_state_array(:,2),aircraft_state_array(:,3),'-.');
hold on;
plot3(aircraft_state_array(1,1),aircraft_state_array(1,2),aircraft_state_array(1,3),"xg");
plot3(aircraft_state_array(1,1),aircraft_state_array(1,2),aircraft_state_array(1,3),"xg");
plot3(aircraft_state_array(end,1),aircraft_state_array(end,2),aircraft_state_array(end,3),"xr")
end



function motor_forces = ComputeMotorForces(Zc,Lc,Mc,Nc,R,km)

c = inv([-1,-1,-1,-1;-R/sqrt(2),-R/sqrt(2),R/sqrt(2),R/sqrt(2); ...
    R/sqrt(2),-R/sqrt(2),-R/sqrt(2),R/sqrt(2); km,-km,km,-km]);

motor_forces = c * [Zc; Lc; Mc; Nc];

end

