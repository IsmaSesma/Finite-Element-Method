% *************************************************************************
%     PROCEDURES TO OBTAIN YOUNG'S MODULUS WITH RESONANCE FREQUENCIES
% *************************************************************************

clc;clear;close all;

%% MASIC AND GEOMETRIC INPUTS

promt.m = 'Mass of the beam (kg): ';
beam.m = input(promt.m);                            % Beam's mass
promt.L = 'Length of the beam (m): ';
beam.L = input(promt.L);                            % Beam's length
promt.b = 'Width of the beam (m): ';
beam.b = input(promt.b);                            % Beam's width
promt.t = 'Thickness of the beam (m): ';
beam.t = input(promt.t);                            % Beam's thickness

beam.rho = beam.m/beam.L/beam.b/beam.t;             % Beam's density
beam.rhom = beam.rho*beam.b*beam.t;                 % Beam's linear mass density
beam.Ixx = beam.b*beam.t^3/12;                      % Beam's area moment of inertia

% λ coefficients for the methods considered

lambda.ISO = [1.87510 4.69410 7.85476];             % ISO-19640 lamdas λ
lambda.AFF = [3.75038 9.39740 15.73438];            % Antiresonances of a Free-Free Beam λ   
lambda.RFF = [4.73004 10.99561 17.27876];           % Resonances of a Free-Free Beam λ 

%% RESONANCE AND ANTI-RESONANCE FREQUENCIES (INPUT FROM DYNAMIC TEST)

promt.rw = 'First 3 resonance frequencies of the beam (Hz): ';
beam.rw = input(promt.rw);
promt.aw = 'First 3 anti-resonance frequencies of the beam (Hz): ';
beam.aw = input(promt.aw);

%% YOUNG'S MODULUS

% ISO - 16940 METHOD
    % The ISO standard works with impedance instead of admitance, so what 
    % they call impedance resonances are what here is considered as
    % admitance antiresonances.

E_ISO = beam.rhom/beam.Ixx*(2*pi.*beam.aw.*(beam.L/2)^2./lambda.ISO.^2).^2*1E-9;
E_ISO_med = mean(E_ISO);                                                            % Mean value
E_ISO_std = std(E_ISO);                                                             % Standard deviation

% ANTIRESONANCES OF A FREE-FREE BEAM METHOD

E_AFF = beam.rhom/beam.Ixx*(2*pi*beam.aw*(beam.L)^2./lambda.AFF.^2).^2*1E-9;
E_AFF_med = mean(E_AFF);                                                            % Mean value
E_AFF_std = std(E_AFF);                                                             % Standard deviation;

% RESONANCES OF A FREE-FREE BEAM METHOD

E_RFF = beam.rhom/beam.Ixx*(2*pi.*beam.rw.*(beam.L)^2./lambda.RFF.^2).^2*1E-9;
E_RFF_med = mean(E_RFF);                                                            % Mean value
E_RFF_std = std(E_RFF);                                                             % Standard deviation

% FREE RESPONSEOF A CANTILEVER BEAM METHOD

disp("Young's Modulus with ISO-16940, Antiresonances of FF Beam and Resonances of FF Beam (GPa)")
disp(num2str([E_ISO', E_AFF', E_RFF']))

