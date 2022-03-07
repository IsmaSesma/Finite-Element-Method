% *************************************************************
%            FEM DYNAMIC ANALYSIS OF 1D FREE-FREE BEAM 
% *************************************************************

clc;clear;close all;
%% MASIC AND GEOMETRIC DATA (all in ISU)

beam.E = 2E9;
beam.L = 0.2;
beam.b = 0.02;
beam.t = 0.004;
beam.Ixx = beam.b*beam.t^3/12;
beam.m = 0.03;
beam.rho = beam.m/beam.L/beam.b/beam.t;
beam.rhoL = beam.m/beam.L;

%% IMPUT DATA

p = [1; 1; 1; 1; 1; 1];         % Input force's amplitudes (each value represents deflection and twist of each node of the beam)
Omega = 1;                      % Frecuency of the input force

%% STIFFNESS AND INERTIA BEAM MATRICES

ne = 2;         % Number of elements to be used (determined by wavelenght and propagation speed of the wave in the beam)
nn = 3;         % Number of nodes
dofn = 2;       % Degrees of freedom per node (only considering flexion)
DOF = dofn*nn;  % Total dof

coord_n = zeros(nn,2);                  % Nodal coordinates matrix
coord_n(:,1) = 0:beam.L/ne:beam.L;      % Only X coordinate is different than 0

connect_e = zeros(ne,2);                % Connectivity matrix of elements through nodes
connect_e(:,1) = 1:1:ne;
connect_e(:,2) = 2:1:nn;

K = zeros(DOF);     % Initilization of the stiffness matrix
M_consist = zeros(DOF); % Consistent mass matrix
M_lumped = zeros(DOF);  % Lumped mass matrix

for e = 1:ne
    index = connect_e(e,:);                                    
    x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
    Le = x2 - x1;                                               % Length of the element
    Je = Le/2; iJe = 1/Je;                                      % Jacobian of the transformation to isoparametric coordinates

    dofe = [index(1)*dofn-1 index(1)*dofn...                    % DOF of each element
            index(2)*dofn-1 index(2)*dofn];

% STIFFNESS
     
    k = beam.E*beam.Ixx/Le^3;
    Kef = k*[12 6*Le -12 6*Le;...
             6*Le 4*Le^2 -6*Le 2*Le^2;...
             -12 -6*Le 12 -6*Le;...
             6*Le 2*Le^2 -6*Le 4*Le^2];

    K(dofe, dofe) = K(dofe, dofe) + Kef;

% INERTIA

    % First compute the consistent mass matrix

    m = beam.rho*beam.b*beam.t*Le;                  % Mass of the element          
    Mce = m/420*[156 22*Le 54 -13*Le; 22*Le 4*Le^2 13*Le -3*Le^2; 54 13*Le 156 -22*Le; -13*Le -3*Le^2 -22*Le 4*Le^2];

    M_consist(dofe,dofe) = M_consist(dofe,dofe) + Mce;

 % Second compute the lumped mass matrix

    Mle = m/2*[1 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0];

    M_lumped(dofe,dofe) = M_lumped(dofe,dofe) + Mle;
end

%% RESOLUTION OF THE DYNAMIC SYSTEM (q0[[K] - w^2[M]] = p0)
% Dumped is assumed negligible

D_c = (K - Omega^2*M_consist);          % Dynamic stiffness matrix with consistent mass matrix
q_c = D_c\p;                            % Displacement's amplitudes with consistent mass matrix

D_l = (K - Omega^2*M_lumped);           % Dynamic stiffness matrix with lumped mass matrix
q_l = D_l\p;                            % Displacement's amplitudes with lumped mass matrix





