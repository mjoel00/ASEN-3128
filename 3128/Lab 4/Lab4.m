%ASEN 3128 Lab IV

%Matlab Maintenance
clear, clc, close all

%Constants
Ix = 5.8*10^-5;     %[kg*m^2] Inertia in the body x
Iy = 7.2*10^-5;     %[kg*m^2] Inertia in the body y
sf = 0.1;           % Scale Factor; relates tau1, tau2
g = 9.81;           %[m/s^2] Acceleration due to Gravity
Tau_1 = 0.5;        %[s] Time constant (given)
Tau_2 = Tau_1*sf;   %[s] Second Time constant (picked by yours truely)

%% Problem 1
%Math
lambda_1 = -1/Tau_1;    % Eigenvalue for tau1
lambda_2 = -1/Tau_2;    % Eigenvalue for tau2

k2_lat = Ix*(lambda_2*lambda_1);
k1_lat = Ix*(-lambda_2-lambda_1);

k2_long = Iy*(lambda_2*lambda_1);
k1_long = Iy*(-lambda_2-lambda_1);

%% Problem 4
% See below for me going on a tangent and being dumb for 30 mins
%{
% Lateral Shenanigans
syms s k1_sym k2_sym k3_sym v_r v g_sym Ix_sym Iy_sym
% Double check my derivation of the transfer function
temp.phi = v*s/g_sym;
temp.p = s*temp.phi;
temp.L_c = Ix_sym*s*temp.p;

%Type up the control law
temp.eq1 = (v_r - v)*k3_sym - temp.phi*k2_sym - temp.p*k1_sym;
temp.eq2 = solve(temp.eq1 == temp.L_c,v);
temp.eq3 = simplify(temp.eq2/v_r);

%Extract the characteristic equation from eq3
%}
% ACTUALLY doing the problem...

%Calculate the constraint for the eigenvalues of the A matricies
T_max = 1.25; %[s] 

%Add a number for the number of K3 values to test
n = 1000;
%-------------------------------------------------------------------------%
figure
subplot(2,1,1)
xline(-1/T_max,'r--');
grid on
hold on %ominous music plays...

%Lateral Shenaningans:
% Store each A matrix et al in the cell thingies below
A_lats = {};
eigvals_lats = {};

% Loops, loops! We love loops! Loop through the K3 possibilities and test
% each of them; storing the resulting A matrix and eigenvalues in the
% respective cell arrays. Also plot :D
%Also... calculate the 
for i = 1:n
    K3_vectorLat(i) = i*Ix/200; % Stolen from prof. Frew
    A_lats{i} = [0,1,0,0;...
                 0,0,g,0;...
                 0,0,0,1;...
                 0,-K3_vectorLat(i)/Ix,-k2_lat/Ix,-k1_lat/Ix];
    %Calculate and store the eigen values for each A matrix
    eigvals_lats{i} = eig(A_lats{i});
    
    %Math time! Find a good K3 Value
    %To start, check if we have all real eigenvalues...
    %... also only do this until we've fount a good value
    if all(imag(eigvals_lats{i}) == 0)
        %Calculate the time constant for each eigenvalue...
        for j = 1:4
            if eigvals_lats{i}(j) == 0
                Ts(j) = 0;
            else
            Ts(j) = -1./eigvals_lats{i}(j);
            end
        end
        %fprintf('Ts: %f,%f,%f,%f\n i: %f\n',Ts(1),Ts(2),Ts(3),Ts(4),i)
        %...Print out the Ts's, by inspection determine the correct index
    end
    %Plotting! Add the four, (x,y)::(R,Im) eigenvalues to the plot
    plot(real(eigvals_lats{i}),imag(eigvals_lats{i}),'bx')
end
%More plotting...
xlim([-2.2,0.1])
title('Locus for Lateral Dynamics')
%By inspection...
k3_latIndex = 376;
K3_lat = K3_vectorLat(k3_latIndex);
%Add the good index to the plot...
plot(real(eigvals_lats{k3_latIndex}),imag(eigvals_lats{k3_latIndex}),'go','LineWidth',1.5)
ylabel('Im')
xlabel('Real')
hold off
%Print it to the command window
fprintf("Lateral K3 Value: %e\n",K3_lat)
%-------------------------------------------------------------------------%
subplot(2,1,2)
xline(-1/T_max,'r--');
grid on
hold on %ominous music plays...

%Longitudinal Shenaningans:
% Store each A matrix et al in the cell thingies below
A_longs = {};
eigvals_longs = {};

% Loops, loops! We love loops! Loop through the K3 possibilities and test
% each of them; storing the resulting A matrix and eigenvalues in the
% respective cell arrays. Also plot :D
%Also... calculate the 
for i = 1:n
    K3_vectorLong(i) = -i*Iy/200; % Stolen from prof. Frew
    A_longs{i} = [0,1,0,0;...
                 0,0,-g,0;...
                 0,0,0,1;...
                 0,-K3_vectorLong(i)/Iy,-k2_long/Iy,-k1_long/Iy];
    %Calculate and store the eigen values for each A matrix
    eigvals_longs{i} = eig(A_longs{i});
    
    %Math time! Find a good K3 Value
    %To start, check if we have all real eigenvalues...
    %... also only do this until we've fount a good value
    if all(imag(eigvals_longs{i}) == 0)
        %Calculate the time constant for each eigenvalue...
        for j = 1:4
            if eigvals_longs{i}(j) == 0
                Ts(j) = 0;
            else
            Ts(j) = -1./eigvals_longs{i}(j);
            end
        end
        %fprintf('Ts: %f,%f,%f,%f\n i: %f\n',Ts(1),Ts(2),Ts(3),Ts(4),i)
        %...Print out the Ts's, by inspection determine the correct index
    end
    %Plotting! Add the four, (x,y)::(R,Im) eigenvalues to the plot
    plot(real(eigvals_longs{i}),imag(eigvals_longs{i}),'bx')
end
%More plotting...
xlim([-2.2,0.1])
title('Locus for Longitudinal Dynamics')
%By inspection...
k3_longIndex = 376;
K3_long = K3_vectorLong(k3_longIndex);
%Add the good index to the plot...
plot(real(eigvals_longs{k3_longIndex}),imag(eigvals_longs{k3_longIndex}),'go','LineWidth',1.5)
ylabel('Im')
xlabel('Real')
hold off
%Print it to the command window
fprintf("Longitudinal K3 Value: %e\n",K3_long)

filename = 'K3_Values';
save(filename,'K3_long','K3_lat');