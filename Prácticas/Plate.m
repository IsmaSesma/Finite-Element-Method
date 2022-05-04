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
dofn = 3;                     % DOF per node
DOF = dofn*nn;                % DOF of the plate

%% NUMERIC INTEGRATION DATA (Gauss-Legendre)

% n = 1 quadrature
chi_ip.f = 0;   % Integration point coordinate
w_ip.f = 2;     % Integration weight
% n = 2 quadrature
chi_ip.s = [-1 1]*1/sqrt(3);   
w_ip.s = [1 1];    
 
%% STIFFNESS AND INERTIA PLATE MATRICES

coord = zeros(nn_s,2);                        % Coordinates of the external nodes
coord(:,1) = 0:plate.a/(nn_s-1):plate.a;      % X coordinates
coord(:,2) = 0:plate.b/(nn_s-1):plate.b;      % Y coordinates

[X,Y] = meshgrid(coord(:,1),coord(:,2));    % Mesh to represent the plate
Z = zeros(nn_s);
surf(X,Y,Z)

K = zeros(nn*dofn);
D = plate.E/(1 - plate.nu)^2*[1 plate.nu 0; plate.nu 1 0; 0 0 (1 - plate.nu)/2];
xe = plate.a; ye = plate.b;
Je = plate.a*plate.b/4;

Ke = zeros(dofn*4);
index = zeros(2);

% Vector with the relative position of every element of the beam
row_col = zeros(ne,2);
for e  = 1:ne
    row_col(e,1) = floor((e-1)/ne_s) + 1;
    row_col(e,2) = mod((e-1),ne_s) + 1;
end

% Connectivity matrix
connectivity = zeros(2,2,ne);
for e = 1:ne
    node1 = (ne_s + 1)*(row_col(e,1) - 1) + row_col(e,2);
    node2 = (ne_s + 1)*(row_col(e,1) - 1) + row_col(e,2) + 1;
    node3 = (ne_s + 1)*row_col(e,1) + row_col(e,2);
    node4 = (ne_s + 1)*row_col(e,1) + row_col(e,2) + 1;

    connectivity(:,:,e) = [node1 node2; node3 node4];

    
end

for e = 1:ne
   index = connectivity(:,:,e);
   x1 = X(index(1,1),1); x2 = X(index(1,2),1);
    for i = 1:2
        for j = 1:2
            chi = chi_ip.s(i); eta = chi_ip.s(j);                   % Value of interpolation points
            w = w_ip.s(i);                                          % Weight of interpolation
            dNdchi_1 = -(1 - eta)/4; dNdeta_1 = -(1 - chi)/4;       % Derivates of form functions
            dNdchi_2 = (1 - eta)/4; dNdeta_2 = -(1 + chi)/4; 
            dNdchi_3 = (1 + eta)/4; dNdeta_3 = (1 + chi)/4; 
            dNdchi_4 = -(1 + eta)/4; dNdeta_4 = (1 - chi)/4; 
            % Kinetic matrix
            B = [0 0 -dNdchi_1*2/xe 0 0 -dNdchi_2*2/xe 0 0 -dNdchi_3*2/xe 0 0 -dNdchi_4*2/xe;
                 0 dNdeta_1*2/ye 0 0 dNdeta_2*2/ye 0 0 dNdeta_3*2/ye 0 0 dNdeta_4*2/ye 0; 
                 0 dNdchi_1*2/xe -dNdeta_1*2/ye 0 dNdchi_2*2/xe -dNdeta_2*2/ye 0 dNdchi_3*2/xe -dNdeta_3*2/ye 0 dNdchi_4*2/xe -dNdeta_4*2/ye]; 
            % Stiffness matrix of the element
            Ke = Ke + B'*D*B*Je*w;
        end
    end
end

%% FUNCTIONS

