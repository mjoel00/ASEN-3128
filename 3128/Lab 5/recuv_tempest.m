% Flight conditions and atmospheric parameters derived from AVL for the 
% University of Colorado's Tempest Unmanned Aircraft
%
%   File created by: Eric Frew, eric.frew@colorado.edu
%   Data taken from files generated by Jason Roadman.
%       - Derivatives come from AVL analysis
%       - Inertias from Solidworks model
%
% If using this data for published work please reference:
%
% Jason Roadman, Jack Elston, Brian Argrow, and Eric W. Frew. 
% “Mission Performance of the Tempest UAS in Supercell Storms.” 
% AIAA Journal of Aircraft, 2012.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% All dimensional parameters in SI units
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aircraft_parameters.g = 9.81;           % Gravitational acceleration [m/s^2]

%Aircraft geometry parameters
aircraft_parameters.S = 0.6282; %[m^2] 
aircraft_parameters.b = 3.067; %[m] 
aircraft_parameters.c = 0.208; %[m]
aircraft_parameters.AR = aircraft_parameters.b*aircraft_parameters.b/aircraft_parameters.S;

% Aircraft mass parameters
aircraft_parameters.m = 5.74; %[kg]
aircraft_parameters.W = aircraft_parameters.m*aircraft_parameters.g; %[N]

% Inertias from Solidworks
SLUGFT2_TO_KGM2 = 14.5939/(3.2804*3.2804);
aircraft_parameters.Ix = SLUGFT2_TO_KGM2*4106/12^2/32.2; %[kg m^2]
aircraft_parameters.Iy = SLUGFT2_TO_KGM2*3186/12^2/32.2; %[kg m^2]
aircraft_parameters.Iz = SLUGFT2_TO_KGM2*7089/12^2/32.2; %[kg m^2]
aircraft_parameters.Ixz = SLUGFT2_TO_KGM2*323.5/12^2/32.2; %[kg m^2]


% Drag terms
aircraft_parameters.e = 0.95; %0.9693; %[-], AVL, Oswald's efficiency factor
aircraft_parameters.K = 1/(pi*aircraft_parameters.AR*aircraft_parameters.e); %drag polar coefficient (CD = CDpa + kCL^2)
aircraft_parameters.CDpa = 0.021; %[-] This is the parasite drag. Total drag is combination of parasite drag and induced drag.
aircraft_parameters.CDmin = 0.021;
aircraft_parameters.CLmin = 0.0;
  
% Engine parameters
aircraft_parameters.Sprop = 0.0707;%0.2027; %[m^2]
aircraft_parameters.Cprop = 1;
aircraft_parameters.kmotor = 40;

% Zero angle of attack aerodynamic forces and moments
aircraft_parameters.CL0 = 0;
aircraft_parameters.Cm0 = 0.1104;

aircraft_parameters.CY0 = 0;
aircraft_parameters.Cl0 = 0;
aircraft_parameters.Cn0 = 0;

% Longtidunal nondimensional stability derivatives from AVL
aircraft_parameters.CLalpha = 6.196683; 
aircraft_parameters.Cmalpha = -1.634010;
aircraft_parameters.CLq = 10.137584; 
aircraft_parameters.Cmq = -24.376066;


% Neglected parameters, check units below if incorporated later
aircraft_parameters.CLalphadot = 0; 
aircraft_parameters.Cmalphadot = 0; 


% Lateral-directional nondimensional stability derivatives from AVL
aircraft_parameters.CYbeta = -0.367231; 
aircraft_parameters.Clbeta = -0.080738; 
aircraft_parameters.Cnbeta = 0.080613; 
aircraft_parameters.CYp = -0.064992; 
aircraft_parameters.Clp = -0.686618; 
aircraft_parameters.Cnp = -0.039384; 
aircraft_parameters.Clr = 0.119718; 
aircraft_parameters.Cnr = -0.052324; 
aircraft_parameters.CYr = 0.213412;  

%Control surface deflection parameters

  % Elevator
  aircraft_parameters.CLde =   0.006776;
  aircraft_parameters.Cmde =  -0.028684; 
  % Aileron
  aircraft_parameters.CYda =  -0.000754;
  aircraft_parameters.Clda =  -0.006290;
  aircraft_parameters.Cnda =  -0.000078;
 
  % Rudder
  aircraft_parameters.CYdr =   0.003056;
  aircraft_parameters.Cldr =   0.000157;
  aircraft_parameters.Cndr =  -0.000856;
  
  
  
  
  
  
