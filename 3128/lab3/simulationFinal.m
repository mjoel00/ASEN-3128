close all
clear
clc

%% Create Constants
global k_p
g = -9.8;
m = 0.068;
r = 0.060;
k_m = 0.0024;
k_p =-0.004;
I_x = 0.000068;
I_y = 0.000092;
I_z = 0.000135;
nu = 0.001;
mu = 0.000002;
col1 = ["ko","b"];              %color arrays used to mark which figures to plot on Non-linear(1) linear(2)
col2 = ["k--","bx"];
col3 = ["k*","b:"];
col4 = ["k","bo"];
col5 = ["kx","b--"];
col6 = ["k:","b*"];
%% Initialize vectors
inertia = [I_x,0,0;0,I_y,0;0,0,I_z];
tspan = [0 5];
deg = 5 * pi/180;
rad = .1;
opts = odeset('Events',@phase);
x_initial = [0;0;10; 0;0;0; 0;0;0; 0;0;0]; % [xe;ye;ze;phi;theta;psi;ue;ve;we;p;q;r]
% Designating state arrays with deviation
xdev1 =[0;0;10; deg;0;0; 0;0;0; 0;0;0];
xdev2 =[0;0;10; 0;deg;0; 0;0;0; 0;0;0];
xdev3 =[0;0;10; 0;0;deg; 0;0;0; 0;0;0];
xdev4 =[0;0;10; 0;0;0; 0;0;0; rad;0;0];
xdev5 =[0;0;10; 0;0;0; 0;0;0; 0;rad;0];
xdev6 =[0;0;10; 0;0;0; 0;0;0; 0;0;rad];
%linear initial systems
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
% [t,state] = ode45( @(t,x) quadFunQ2(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor), tspan, x_initial);
%ode 45 models with each deviation
[t1,state1] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor), tspan, xdev1,opts);
[t2,state2] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor), tspan, xdev2,opts);
[t3,state3] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor), tspan, xdev3,opts);
[t4,state4] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor), tspan, xdev4,opts);
[t5,state5] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor), tspan, xdev5,opts);
[t6,state6] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor), tspan, xdev6,opts);
%% Plotting first deviations column 1 of col designating which figure to plot on li
PlotAircraftSim(t1,state1,col1(1));
PlotAircraftSim(t2,state2,col2(1));
PlotAircraftSim(t3,state3,col3(1));
PlotAircraftSim(t4,state4,col4(1));
PlotAircraftSim(t5,state5,col5(1));
PlotAircraftSim(t6,state6,col6(1));
F_c =  [x_initial(10),x_initial(11),x_initial(12)].*k_p;

PlotAircraftSim(j,y,col1(2))



%% Problem 5
col7 = ["m:","c:"];
col8 = ["m--","co"];
col9 = ["m.-","c.-"];
xdev5_4 =[0;0;10; 0;0;0; 0;0;0; rad;0;0];
xdev5_5 =[0;0;10; 0;0;0; 0;0;0; 0;rad;0];
xdev5_6 =[0;0;10; 0;0;0; 0;0;0; 0;0;rad];

[t5_4,state5_4] = ode45( @(t,x) quadFunQ2(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor), tspan, xdev5_4,opts);
[t5_5,state5_5] = ode45( @(t,x) quadFunQ2(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor), tspan, xdev5_5,opts);
[t5_6,state5_6] = ode45( @(t,x) quadFunQ2(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor), tspan, xdev5_6,opts);

% PlotAircraftSim(t5_4,state5_4,col7(1));
% PlotAircraftSim(t5_5,state5_5,col8(1));
% PlotAircraftSim(t5_6,state5_6,col9(1));

% for a = 1:36
%     figure(a)
%     legend('Non-Linear','Linear)
% end
% figure(5)
% legend('Roll 5 deg','','','Pitch 5 deg','','','Yaw 5 deg','','','Roll Rate .1','','',...
%     'Pitch Rate .1','','','Yaw Rate .1','','','Sim','','','kp Roll','','','kp Pitch','','','kp Yaw')

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
% array location denotes which figure to plot on
% subplots
%creating arrays of labels for use
titles = ["X Postion","Y Postion","Z Postion","Euler Angle Phi","Euler Angle Theta","Euler Angle Psi",...
    "X Velocity","Y Velocity","Z Velocity","Angular Velocity P","Angular Velocity Q","Angular Velocity R"];
xlabels = "Time";
ylabels = [("[m]"),("[m]"),("[m]"),("[rad]"),("[rad]"),("[rad]"),("[m/s]"),...
    ("[m/s]"),("[m/s]"),("[rad/s]"),("[rad/s]"),("[rad/s]")];

% col1 = ["ko","b"];
% col2 = ["k--","bx"];
% col3 = ["k*","b:"];
% col4 = ["k","bo"];
% col5 = ["kx","b--"];
% col6 = ["k:","b*"];

% col7 = ["m:","c:"];
% col8 = ["m--","co"];
% col9 = ["m.-","c.-"]
% array will change figure starting location based on col code found in
% input
if strcmp(col,'ko')||strcmp(col,'b')
    n = 1;
elseif strcmp(col,'k--')||strcmp(col,'bx')
    n = 7;
elseif strcmp(col,'k*')||strcmp(col,'b:')
    n = 13;
elseif strcmp(col,'k')||strcmp(col,'b:')
    n = 19;
elseif strcmp(col,'kx')||strcmp(col,'b--')
    n = 25;
elseif strcmp(col,'k:')||strcmp(col,'b*')
    n = 31;
elseif strcmp(col,'m:')||strcmp(col,'c:')
    n = 37;
elseif strcmp(col,'m--')||strcmp(col,'c--')
    n = 43;
else
    n = 49;
end
figure(n);
subplot(3,1,1);
hold on;
plot(time,aircraft_state_array(:,1),col);
%             %xlim([0 5])
xlabel(xlabels);
ylabel(ylabels(1));
title(titles(1));

subplot(3,1,2);
hold on;
plot(time,aircraft_state_array(:,2),col);
%             %xlim([0 5])
xlabel(xlabels);
ylabel(ylabels{2});
title(titles(2));


subplot(3,1,3);
hold on;
plot(time,aircraft_state_array(:,3),col);
%             %xlim([0 5])
xlabel(xlabels);
ylabel(ylabels{3});
title(titles(3));

figure(n+1);
subplot(3,1,1);
hold on;
plot(time,aircraft_state_array(:,4),col);
%xlim([0 5])
xlabel(xlabels);
ylabel(ylabels(4));
title(titles(4));

subplot(3,1,2);
hold on;
plot(time,aircraft_state_array(:,5),col);
%xlim([0 5])
xlabel(xlabels);
ylabel(ylabels{5});
title(titles(5));


subplot(3,1,3);
hold on;
plot(time,aircraft_state_array(:,6),col);
%xlim([0 5])
xlabel(xlabels);
ylabel(ylabels{6});
title(titles(6));

figure(n+2);
subplot(3,1,1);
hold on;
plot(time,aircraft_state_array(:,7),col);
%xlim([0 5]
xlabel(xlabels);
ylabel(ylabels(7));
title(titles(7));

subplot(3,1,2);
hold on;
plot(time,aircraft_state_array(:,8),col);
%xlim([0 5])
xlabel(xlabels);
ylabel(ylabels{8});
title(titles(8));


subplot(3,1,3);
hold on;
plot(time,aircraft_state_array(:,9),col);
%xlim([0 5])
xlabel(xlabels);
ylabel(ylabels{9});
title(titles(9));

figure(n+3);
subplot(3,1,1);
hold on;
plot(time,aircraft_state_array(:,10),col);
%xlim([0 5])
xlabel(xlabels);
ylabel(ylabels(10));
title(titles(10));

subplot(3,1,2);
hold on;
plot(time,aircraft_state_array(:,11),col);
%xlim([0 5])
xlabel(xlabels);
ylabel(ylabels{11});
title(titles(11));


subplot(3,1,3);
hold on;
plot(time,aircraft_state_array(:,12),col);
%xlim([0 5])
xlabel(xlabels);
ylabel(ylabels{12});
title(titles(12));

figure(n+4)
plot3(aircraft_state_array(:,1),aircraft_state_array(:,2),aircraft_state_array(:,3));
hold on;
plot3(aircraft_state_array(1,1),aircraft_state_array(1,2),aircraft_state_array(1,3),'go','markersize',10);
plot3(aircraft_state_array(end,1),aircraft_state_array(end,2),aircraft_state_array(end,3),'xr','markersize',10)
xlabel('X distance');
ylabel('Y distance');
zlabel('Z distance');
title('Flight path');

        %% input array left commented out until we can generate the array
%             figure(n+5)
%             %xlim([0 5])
%             subplot(4,1,3);
%             plot(time,control_input_array(:,1));
%             hold on;
%             xlabel(xlabels);
%             ylabel("Input 1");
%             title("Input 1");
%             
%             subplot(4,1,2);
%             %xlim([0 5])
%             plot(time,control_input_array(:,2));
%             hold on;
%             xlabel(xlabels);
%             ylabel("Input 2");
%             title("Input 2");
%             
%             subplot(3,1,3);
%             %xlim([0 5])
%             plot(time,control_input_array(:,3));
%             hold on;
%             xlabel(xlabels);
%             ylabel("Input 3");
%             title("Input 3");
%             
%             subplot(3,1,4);
%             %xlim([0 5])
%             plot(time,control_input_array(:,4));
%             hold on;
%             xlabel(xlabels);
%             ylabel("Input 4");
%             title("Input 4");

            
          
    
    
end







function motor_forces = ComputeMotorForces(Zc,Lc,Mc,Nc,R,km)

c = ([-1,-1,-1,-1;-R/sqrt(2),-R/sqrt(2),R/sqrt(2),R/sqrt(2); ...
    R/sqrt(2),-R/sqrt(2),-R/sqrt(2),R/sqrt(2); km,-km,km,-km]);

motor_forces = c\[Zc; Lc; Mc; Nc];

end

function [value, istermainal, direction] = phase(~,x)
    
value = 0>=x(3);
istermainal = 1;
direction = 1;
end