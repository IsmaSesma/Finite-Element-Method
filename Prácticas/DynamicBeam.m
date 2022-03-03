% *************************************************************
%            FEM DYNAMIC ANALYSIS OF 1D FREE-FREE BEAM 
% *************************************************************
clc;clear;close all;
%% MASIC AND GEOMETRIC DATA (all in ISU)
beam.E = 20E9;
beam.L = 1;
beam.b = 1;
beam.t = 1;
beam.Ixx = beam.b*beam.t^3/12;
beam.m = 1;
beam.rho = beam.m/beam.L/beam.b/beam.t;
beam.rhoL = beam.m/beam.L;

%% NUMERIC INTEGRATION DATA (used for K and M matrices)
