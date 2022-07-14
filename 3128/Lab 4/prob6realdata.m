%% Problem 6 Data Analytics
%Gives us the motor rates, time vector, and state vector of the actual
%quadrotor

clc; clear all; close all;

%Load in the data file (must be in the same folder as this script!)
rawData = load('RSdata_10_22_1401_5.mat');
%Time 
REAL_t = rawData.rt_estim.time(:);
%Extract motor command rates
for i = 1:length(rawData.rt_motor.signals.values) %For each measurement... (row)
    for j = 1:4 %For each motor... (column)
        if rawData.rt_motor.signals.values(i,j) < 0 %if the motor force is negative...
            REAL_motorRates(i,j) = -1*sqrt(abs(rawData.rt_motor.signals.values(i,j))*13840.4);
        else %if the motor force is positive...
            REAL_motorRates(i,j) = sqrt(rawData.rt_motor.signals.values(i,j)*13840.4);
        end
    end
end
%Extract state vector
REAL_X = rawData.rt_estim.signals.values;
%% Plotting these...
fig_a = [1,2,3,4,5,6]; % figure a number vector
PlotAircraftSim(REAL_t,REAL_X,'Real Data','b-',fig_a,REAL_motorRates)