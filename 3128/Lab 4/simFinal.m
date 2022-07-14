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
root = sqrt(2);
F = [-1,-1,-1,-1;...
    (-r/root),(-r/root),r/root,r/root;...
     r/root,(-r/root),(-r/root),r/root;...
     k_m,(-k_m),k_m,(-k_m)]^-1 * [m*g;0;0;0];
f1 = F(1); f2 = F(2); f3 = F(3); f4 = F(4);
F = [f1;f2;f3;f4];
col1 = ["ko","b"];              %color arrays used to mark which figures to plot on Non-linear(1) linear(2)
col2 = ["ks","bx"];
col3 = ["k*","bs"];
col4 = ["k","bo"];
col5 = ["kx","bd"];
col6 = ["kd","b*"];
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
[y1,j1,x1] = initial(sys,xdev1);
%% Computing motor forces
mf31 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf31 = (zeros(length(j1),4)+1).*mf31;
[y2,j2,x2] = initial(sys,xdev2);
mf32 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf32 = (zeros(length(j2),4)+1).*mf32;
[y3,j3,x3] = initial(sys,xdev3);
mf33 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf33 = (zeros(length(j3),4)+1).*mf33;
[y4,j4,x4] = initial(sys,xdev4);
mf34 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf34 = (zeros(length(j4),4)+1).*mf34;
[y5,j5,x5] = initial(sys,xdev5);
mf35 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf35 = (zeros(length(j5),4)+1).*mf35;
[y6,j6,x6] = initial(sys,xdev6);
mf36 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf36 = (zeros(length(j6),4)+1).*mf36;

%% setting motor forces for start state
L_c = x_initial(10) *  k_p;
M_c = x_initial(11) * k_p;
N_c = x_initial(12) * k_p;
forces_motor = ComputeMotorForces(m*g,L_c,M_c,N_c,r,k_m);

%% ode 45 models with each deviation then computing motor forces for duration
[t1,state1] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor'), tspan, xdev1,opts);
mf1 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf1 = (zeros(length(t1),4)+1).*mf1; % creating force array 

[t2,state2] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor'), tspan, xdev2,opts);
mf2 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf2 = (zeros(length(t2),4)+1).*mf2;

[t3,state3] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor'), tspan, xdev3,opts);
mf3 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf3 = (zeros(length(t3),4)+1).*mf3;

[t4,state4] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor'), tspan, xdev4,opts);
mf4 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf4 = (zeros(length(t4),4)+1).*mf4;

[t5,state5] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor'), tspan, xdev5,opts);
mf5 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf5 = (zeros(length(t5),4)+1).*mf5;

[t6,state6] = ode45( @(t,x) quadFunQ1(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor'), tspan, xdev6,opts);
mf6 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf6 = (zeros(length(t6),4)+1).*mf6;

%% Plotting first deviations column 1 of col designating which figure to plot on
PlotAircraftSim(t1,state1,mf1,col1(1));
PlotAircraftSim(t2,state2,mf2,col2(1));
PlotAircraftSim(t3,state3,mf3,col3(1));
PlotAircraftSim(t4,state4,mf4,col4(1));
PlotAircraftSim(t5,state5,mf5,col5(1));
PlotAircraftSim(t6,state6,mf6,col6(1));
F_c =  [x_initial(10),x_initial(11),x_initial(12)].*k_p;

%% including the same motor forces since neither linear or non linear models are changing motor forces
PlotAircraftSim(j1,y1,mf31,col1(2));
PlotAircraftSim(j2,y2,mf32,col2(2));
PlotAircraftSim(j3,y3,mf33,col3(2));
PlotAircraftSim(j4,y4,mf34,col4(2));
PlotAircraftSim(j5,y5,mf35,col5(2));
PlotAircraftSim(j6,y6,mf36,col6(2));



%% Problem 5 New col for designation and id of kp plots
col7 = ["m*","c:"];
col8 = ["ms","co"];
col9 = ["mx","c.-"];
xdev5_4 =[0;0;10; 0;0;0; 0;0;0; rad;0;0];
xdev5_5 =[0;0;10; 0;0;0; 0;0;0; 0;rad;0];
xdev5_6 =[0;0;10; 0;0;0; 0;0;0; 0;0;rad];

[t5_4,state5_4] = ode45( @(t,x) quadFunQ2(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor'), tspan, xdev5_4,opts);
Zc54 = (m*g) * ones(length(t5_4),1);
Lc54 = -k_m*state5_4(:,10)-k_m*state5_4(:,4);  % creating control moments based of output state array
Mc54 = -k_m*state5_4(:,11)-k_m*state5_4(:,5);
Nc54 = -k_m*state5_4(:,12)-k_m*state5_4(:,6);
mf54 = ComputeMotorForces(m*g,Lc54,Mc54,Nc54,r,k_m);
mf54 = (zeros(length(t5_4),4)+1).*mf54;

[t5_5,state5_5] = ode45( @(t,x) quadFunQ2(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor'), tspan, xdev5_5,opts);
Zc54 = (m*g) * ones(length(t5_5),1);
Lc55 = -k_m*state5_5(:,10)-k_m*state5_5(:,4);
Mc55 = -k_m*state5_5(:,11)-k_m*state5_5(:,5);
Nc55 = -k_m*state5_5(:,12)-k_m*state5_5(:,6);
mf55 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf55 = (zeros(length(t5_5),4)+1).*mf55;

[t5_6,state5_6] = ode45( @(t,x) quadFunQ2(t,x,g,m,inertia,k_m,nu,mu,r,forces_motor'), tspan, xdev5_6,opts);
Zc56 = (m*g) * ones(length(t5_6),1);
Lc56 = -k_m*state5_6(:,10)-k_m*state5_6(:,4);
Mc56 = -k_m*state5_6(:,11)-k_m*state5_6(:,5);
Nc56 = -k_m*state5_6(:,12)-k_m*state5_6(:,6);
mf56 = ComputeMotorForces(m*g,0,0,0,r,k_m);
mf56 = (zeros(length(t5_6),4)+1).*mf56;

%% adding plots to figures for comparison
PlotAircraftSim(t5_4,state5_4,mf54,col7(1));
PlotAircraftSim(t5_5,state5_5,mf55,col8(1));
PlotAircraftSim(t5_6,state5_6,mf56,col9(1));
%% setting main titles for plots
names = ["5 Deg Roll","5 Deg Pitch","5 Deg Yaw",".1 Roll Rate",".1 Pitch Rate",".1 Yaw Rate"];
loop = 1;
for a = 1:6:36
    if a <=18
    figure(a)
    legend('Non-Linear','Linear')
    sgtitle(names(loop));
    figure(a+1)
    legend('Non-Linear','Linear')
    sgtitle(names(loop));
    figure(a+2)
    legend('Non-Linear','Linear')
    sgtitle(names(loop));
    figure(a+3)
    legend('Non-Linear','Linear')
    sgtitle(names(loop));
    figure(a+4)
    legend('Non-Linear','','','Linear','','')
    sgtitle(names(loop));
    loop=loop+1;
    figure(a+5)
    sgtitle(names(loop));
    legend('Non-Linear','Linear')

    else
    figure(a)
    legend('Non-Linear','Linear','Kp')
    sgtitle(names(loop));
    figure(a+1)
    legend('Non-Linear','Linear','Kp')
    sgtitle(names(loop));
    figure(a+2)
    legend('Non-Linear','Linear','Kp')
    sgtitle(names(loop));
    figure(a+3)
    legend('Non-Linear','Linear','Kp')
    sgtitle(names(loop));
    figure(a+4)
    legend('Non-Linear','','','Linear','','','Kp','','')
    sgtitle(names(loop));
    
    figure(a+5)
    legend('Non-Linear','Linear','Kp')
    sgtitle(names(loop));
    loop = loop+1;

        
    end

        
end

%% autosaving plots since publising cuts off the titles toward the latter plots.
% % a pop up will appear to save plots to location, enter name and location.
% n=36;
% [name,path] = uiputfile('*.*');
% filename = [path,name];
% plotname = cell(n,1);
% for i = 1:36
%     plotname{i} = [filename,'_',sprintf('%02.0f',i),'.png'];
%     fig_handle = figure(i);
%     saveas(fig_handle,plotname{i},'png');
% end


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


function [x,t] = PlotAircraftSim(time,aircraft_state_array,control_input_array,col)
% array location denotes which figure to plot on
% subplots
%creating arrays of labels for use
titles = ["X Postion","Y Postion","Z Postion","Euler Angle Phi","Euler Angle Theta","Euler Angle Psi",...
    "X Velocity","Y Velocity","Z Velocity","Angular Velocity P","Angular Velocity Q","Angular Velocity R"];
xlabels = "Time";
ylabels = [("[m]"),("[m]"),("[m]"),("[rad]"),("[rad]"),("[rad]"),("[m/s]"),...
    ("[m/s]"),("[m/s]"),("[rad/s]"),("[rad/s]"),("[rad/s]")];

% col1 = ["ko","b"];   %the cols used as reference to build the strcmp statements
% col2 = ["ks","bx"];
% col3 = ["k*","bs"];
% col4 = ["k","bo"];
% col5 = ["kx","bd"];
% col6 = ["kd","b*"];
% col7 = ["m*","c:"];
% col8 = ["ms","co"];
% col9 = ["mx","c.-"];
% array will change figure starting location based on col code found in
% input
if strcmp(col,'ko')||strcmp(col,'b')
    n = 1;
elseif strcmp(col,'ks')||strcmp(col,'bx')
    n = 7;
elseif strcmp(col,'k*')||strcmp(col,'bs')
    n = 13;
elseif strcmp(col,'k')||strcmp(col,'bo')
    n = 19;
elseif strcmp(col,'kx')||strcmp(col,'bd')
    n = 25;
elseif strcmp(col,'kd')||strcmp(col,'b*')
    n = 31;
elseif strcmp(col,'m*')||strcmp(col,'c:')
    n = 19;
elseif strcmp(col,'ms')||strcmp(col,'c--')
    n = 25;
else
    n = 31;
end
figure(n);
subplot(3,1,1);
hold on;
plot(time,aircraft_state_array(:,1),col);
xlim([0 2])
xlabel(xlabels);
ylabel(ylabels(1));
title(titles(1));

subplot(3,1,2);
hold on;
plot(time,aircraft_state_array(:,2),col);
xlim([0 2])
xlabel(xlabels);
ylabel(ylabels{2});
title(titles(2));


subplot(3,1,3);
hold on;
plot(time,aircraft_state_array(:,3),col);
xlim([0 2])
xlabel(xlabels);
ylabel(ylabels{3});
title(titles(3));

figure(n+1);
subplot(3,1,1);
hold on;
plot(time,aircraft_state_array(:,4),col);
xlim([0 2])
xlabel(xlabels);
ylabel(ylabels(4));
title(titles(4));

subplot(3,1,2);
hold on;
plot(time,aircraft_state_array(:,5),col);
xlim([0 2])
xlabel(xlabels);
ylabel(ylabels{5});
title(titles(5));


subplot(3,1,3);
hold on;
plot(time,aircraft_state_array(:,6),col);
xlim([0 2])
xlabel(xlabels);
ylabel(ylabels{6});
title(titles(6));

figure(n+2);
subplot(3,1,1);
hold on;
plot(time,aircraft_state_array(:,7),col);
xlim([0 2])
xlabel(xlabels);
ylabel(ylabels(7));
title(titles(7));

subplot(3,1,2);
hold on;
plot(time,aircraft_state_array(:,8),col);
xlim([0 2])
xlabel(xlabels);
ylabel(ylabels{8});
title(titles(8));


subplot(3,1,3);
hold on;
plot(time,aircraft_state_array(:,9),col);
xlim([0 2])
xlabel(xlabels);
ylabel(ylabels{9});
title(titles(9));

figure(n+3);
subplot(3,1,1);
hold on;
plot(time,aircraft_state_array(:,10),col);
xlim([0 2])
xlabel(xlabels);
ylabel(ylabels(10));
title(titles(10));

subplot(3,1,2);
hold on;
plot(time,aircraft_state_array(:,11),col);
xlim([0 2])
xlabel(xlabels);
ylabel(ylabels{11});
title(titles(11));


subplot(3,1,3);
hold on;
plot(time,aircraft_state_array(:,12),col);
xlim([0 2])
xlabel(xlabels);
ylabel(ylabels{12});
title(titles(12));

figure(n+4)
plot3(aircraft_state_array(:,1),aircraft_state_array(:,2),aircraft_state_array(:,3),col);
hold on;
plot3(aircraft_state_array(1,1),aircraft_state_array(1,2),aircraft_state_array(1,3),'go','markersize',10);
plot3(aircraft_state_array(end,1),aircraft_state_array(end,2),aircraft_state_array(end,3),'xr','markersize',10)
xlabel('X distance');
ylabel('Y distance');
zlabel('Z distance');
title('Flight path');

        %% input array left commented out until we can generate the array
figure(n+5)
subplot(4,1,1)
plot(time,control_input_array(:,1),col);
xlim([0 2])
hold on;
xlabel(xlabels);
ylabel("[N]");
title("Input 1");

subplot(4,1,2);
xlim([0 2])
plot(time,control_input_array(:,2),col);
hold on;
xlabel(xlabels);
ylabel("[N]");
title("Input 2");

subplot(4,1,3);
xlim([0 2])
plot(time,control_input_array(:,3),col);
hold on;
xlabel(xlabels);
ylabel("[N]");
title("Input 3");

subplot(4,1,4);
xlim([0 2])
plot(time,control_input_array(:,4),col);
hold on;
xlabel(xlabels);
ylabel("[N]");
title("Input 4");

            
          
    
    
end



function motor_forces = ComputeMotorForces(Zc,Lc,Mc,Nc,R,km)
motor_forces = zeros(length(Zc),4);
c = [-1,-1,-1,-1;...
    -R/sqrt(2),-R/sqrt(2),R/sqrt(2),R/sqrt(2); ...
    R/sqrt(2),-R/sqrt(2),-R/sqrt(2),R/sqrt(2);...
    km,-km,km,-km]^-1;
for i = 1:length(Zc)
    motor_forces(i,:) = c * [Zc(i);Lc(i);Mc(i);Nc(i)];

end

end

function [value, istermainal, direction] = phase(~,x)
    
value = 0>=x(3);
istermainal = 1;
direction = 1;
end