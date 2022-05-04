% ***********************************************************************
%           FEM MODEL OF A PLATE USING RECTANGULAR ELEMENTS 
% ***********************************************************************

clc;clear;close all;

%% MASIC AND GEOMETRIC DATA

plate.E = 2E9; plate.nu = 0.3;
plate.a = 0.2; plate.b = 0.2;                % Plate's dimensions
plate.t = 0.04;                              % Plate's thickness
plate.m = 0.5;                               % Plate's mass

%% INPUT DATA

ne_s = 3;                     % Number of bidimensional elements per side of the plate
nn_s = ne_s+1;                % Number of nodes per side
ne = ne_s^2;                  % Number of bidimensional elements in the whole plate
nn = nn_s^2;                  % Number of nodes
dofn = 2;                     % DOF per node
DOF = dofn*nn;                % DOF of the plate

%% STIFFNESS AND INERTIA PLATE MATRICES

coord_n = zeros(nn_s,2);                        % Nodal coordinates matrix
coord_n(:,1) = 0:plate.a/(nn_s-1):plate.a;      % X coordinates
coord_n(:,2) = 0:plate.b/(nn_s-1):plate.b;      % Y coordinates

[X,Y] = meshgrid(coord_n(:,1),coord_n(:,2));    % Mesh to represent the plate
Z = zeros(nn_s);
surf(X,Y,Z)

connect = 1:1:nn;                             % Connectivity matrix of the elements through the nodes
connect_e = reshape(connect,[4,4])';

K = zeros(nn*dofn);

for e = 1:ne
    index = 
end

%% FUNCTIONS

function N_f = N(chi, eta)              % Form functions of the bidimensional element of 4 nodes
    N_f = 1/4*(1-chi)*(1-eta);
end