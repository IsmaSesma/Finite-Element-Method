% *************************************************************
%            FEM DYNAMIC ANALYSIS OF 1D FREE-FREE BEAM 
% *************************************************************

% This file works as a compendium where I have been coding
% different methods and alternatives to do the same things, 
% for the final and optimize code go DynamicBeam_def.

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

%% NUMERIC INTEGRATION DATA (Gauss-Legendre)

% n = 1 quadrature
chi_ip.f = 0;   % Integration point coordinate
w_ip.f = 2;     % Integration weight
% n = 2 quadrature
chi_ip.s = [-1 1]*1/sqrt(3);   
w_ip.s = [1 1];     
% n = 3 quadrature
chi_ip.t = [-sqrt(3/5) 0 sqrt(3/5)];   
w_ip.t = [5/9 8/9 5/9];   
% n = 4 quadrature
chi_ip.ft = [-sqrt((3+2*sqrt(6/5))/7) -sqrt((3-2*sqrt(6/5))/7)  sqrt((3-2*sqrt(6/5))/7) sqrt((3+2*sqrt(6/5))/7)];   
w_ip.ft = [(18-sqrt(30))/36 (18-sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36];
% n = 5 quadrature
chi_ip.fv = [-1/3*sqrt(5+2*sqrt(10/7))  -1/3*sqrt(5-2*sqrt(10/7)) 0 1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7))];   
w_ip.fv = [(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900]; 

%% STIFFNESS BEAM MATRIX

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
K_ = zeros(DOF);

for e = 1:ne
    index = connect_e(e,:);                                    
    x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
    Le = x2 - x1;                                               % Length of the element
    Je = Le/2; iJe = 1/Je;                                      % Jacobian of the transformation to isoparametric coordinates

    dofe = [index(1)*dofn-1 index(1)*dofn...                    % DOF of each element
            index(2)*dofn-1 index(2)*dofn];

    Kef = zeros(4);     % Initialization of the flexion stiffness matrix of one element
    for i = 1:2
        chi = chi_ip.s(i);
        w = w_ip.s(i);
        N = sd_Nh(chi,Je);       % Second derivative of hermitian form functions
        B = iJe^2*N;             % Bending kinetic matrix

        Kef = Kef + B'*beam.E*beam.Ixx*B*Je*w;
    end

    K(dofe, dofe) = K(dofe, dofe) + Kef;

     % Since the beam is a simple structrue to model, it can also be done with
     % the analytical formula available in literature

     Kef_ = zeros(4);
     
     k = beam.E*beam.Ixx/Le^3;
     Kef_ = k*[12 6*Le -12 6*Le;...
               6*Le 4*Le^2 -6*Le 2*Le^2;...
               -12 -6*Le 12 -6*Le;...
               6*Le 2*Le^2 -6*Le 4*Le^2];

     K_(dofe, dofe) = K_(dofe, dofe) + Kef_;
end

%% INERTIA BEAM MATRIX

M_consist = zeros(DOF); % Consistent mass matrix
M_lumped = zeros(DOF);  % Lumped mass matrix
M_consist_2 = zeros(DOF);

for e = 1:ne
    index = connect_e(e,:);                                    
    x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
    Le = x2 - x1;                                              
    Je = Le/2; iJe = 1/Je;                                      

    dofe = [index(1)*dofn-1 index(1)*dofn...                   
            index(2)*dofn-1 index(2)*dofn];

% First compute the consistent mass matrix

    Mce = zeros(4);     % Initialization of the consistent mass matrix of one element
    for i = 1:5
        chi = chi_ip.fv(i);
        w = w_ip.fv(i);
        N = Nh(chi,Je);       % Hermitian form functions

        Mce = Mce + N'*beam.rho*beam.b*beam.t*N*Je*w;
    end

    M_consist(dofe,dofe) = M_consist(dofe,dofe) + Mce;

 % Since the beam is a simple structrue to model, it can also be done with
 % the analytical formula available in literature

    k = beam.rho*beam.b*beam.t*Le/420;
    Mce_ = k*[156 22*Le 54 -13*Le; 22*Le 4*Le^2 13*Le -3*Le^2; 54 13*Le 156 -22*Le; -13*Le -3*Le^2 -22*Le 4*Le^2];

    M_consist_2(dofe,dofe) = M_consist_2(dofe,dofe) + Mce_;         % Will be using this one to simplify computation

 % Second compute the lumped mass matrix

    Mle = zeros(4);     % Initialization of the lumped mass matrix of one element
    Mle = beam.rho*beam.b*beam.t*1/nn*[1/3 0 0 0; 0 0 0 0; 0 0 1/3 0; 0 0 0 0];

    M_lumped(dofe,dofe) = M_lumped(dofe,dofe) + Mle;
end


















%% FUNCTIONS

function NH = Nh(chi,Je)          % Hermitic form functions
    NH = 1/4*[(2-3*chi+chi^3) Je*(1-chi-chi^2+chi^3) (2+3*chi-chi^3) Je*(-1-chi+chi^2+chi^3)];
end

function sd_NH = sd_Nh(chi,Je)    % Second derivative of hermitic form functions
    sd_NH = [(chi*3/2) Je*(chi*3-1)/2 -(chi*3/2) Je*(chi*3+1)/2];
end

