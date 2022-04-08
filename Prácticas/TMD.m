% ***********************************************************************
%           ANALYSIS OF A CANTILEVER BEAM AS A TUNNED MASS DUMPER 
% ***********************************************************************

clc;clear;close all;

%% MASIC AND GEOMETRIC DATA (all in ISU)

beam.E = 2E9;                                       % Elastic modulus
beam.L = 0.2;                                       % Beam's length
beam.b = 0.02;                                      % Beam's width
beam.t = 0.004;                                     % Beam's thickness
beam.tl = 0.02;                                     % Thickness of the last element of the beam      
beam.m = 0.03;                                      % Beam's mass
beam.Ixx = beam.b*beam.t^3/12;                      % Beam's area moment of inertia
beam.rho = beam.m/beam.L/beam.b/beam.t;             % Beam's density
beam.rhom = beam.rho*beam.b*beam.t;                 % Beam's linear mass density
beam.x = (0:0.001:beam.L);                          % Beam's partition 
a = 1E-6; b = 1E-6;                                 % Propotional damping coefficient

%% INPUT DATA

ne = 100;                       % Number of elements to be used (determined by wavelenght and propagation speed of the wave in the beam)
nn = ne + 1;                    % Number of nodes
dofn = 2;                       % Degrees of freedom per node (only considering flexion)
DOF = dofn*nn;                  % Total dof
RDOF = 2;                       % Number of restricted DOF
FDOF = DOF - RDOF;              % Number of free DOF

p = zeros(DOF,1);
p(DOF/2)= 1;                    % Input force's amplitudes (each value represents deflection and twist of each node of the beam)
F = 2000;                       % Maximum frecuency
f = (1:1:F);                    % Frecuency sweep of the input force

beam.modes = (1:7);             % Shape mode plotted (ascending order: first 2 are rigid solid modes and the other are the ones we want)
beam.resonance_i = [3 5 7];     % Index i is where the resonance frequencies are according to the theory

free_dof = (3:DOF);             % Free DOF                    
restr_dof = [1 2];              % Restricted DOF in the fixed end

%% STIFFNESS AND INERTIA BEAM MATRICES

coord_n = zeros(nn,2);                  % Nodal coordinates matrix
coord_n(:,1) = 0:beam.L/ne:beam.L;      % Only X coordinate is different than 0

connect_e = zeros(ne,2);                % Connectivity matrix of elements through nodes
connect_e(:,1) = 1:1:ne;
connect_e(:,2) = 2:1:nn;

K = zeros(DOF);                         % Initilization of the stiffness matrix
M = zeros(DOF);                         % Consistent mass matrix

for e = 1:ne
    index = connect_e(e,:);                                    
    x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
    Le = x2 - x1;                                               % Length of the element

    dofe = [index(1)*dofn-1 index(1)*dofn...                    % DOF of each element
            index(2)*dofn-1 index(2)*dofn];
            
%     if e == ne                                                  % Changes the thickness of the beam in order to simulate a TMD
%         beam.t = beam.tl;
%         beam.Ixx = beam.b*beam.t^3/12;
%         beam.rho = beam.m/beam.L/beam.b/beam.t;
%     else
%         beam.t = beam.t;
%     end

% STIFFNESS
     
    k = beam.E*beam.Ixx/Le^3;
    Kef = k*[12 6*Le -12 6*Le;...
             6*Le 4*Le^2 -6*Le 2*Le^2;...
             -12 -6*Le 12 -6*Le;...
             6*Le 2*Le^2 -6*Le 4*Le^2];

    K(dofe, dofe) = K(dofe, dofe) + Kef;

% INERTIA

  % Consistent mass matrix

    m_e = beam.rho*beam.b*beam.t*Le;            % Mass of the element          
    Mce = m_e/420*[156 22*Le 54 -13*Le;...      % Consistent mass matrix of each element
          22*Le 4*Le^2 13*Le -3*Le^2;...
          54 13*Le 156 -22*Le;...
          -13*Le -3*Le^2 -22*Le 4*Le^2];     

    M(dofe,dofe) = M(dofe,dofe) + Mce;

end

K_FF = K(free_dof, free_dof);                   % Stiffness matrix of the free-free DOF
M_FF = M(free_dof, free_dof);                   % Inertia matrix of the free-free DOF
p_F = p(free_dof);                              % Load's vector in the free DOF

%% RESOLUTION OF THE CONSERVATIVE DYNAMIC SYSTEM (q0[[K] - Ω^2[M]] = p0)

% Dumped is assumed negligible

D_FF = zeros(FDOF,FDOF,F);
q0_c = zeros(FDOF,F);

for i = 1:F           % This loop makes a sweep in the frecuencies up to the maximum frecuency of interest (Hz)

    D_FF(:,:,i) = (K_FF - (2*pi*i)^2*M_FF);                                          % Dynamic stiffness matrix with consistent mass matrix
    q0_c(:,i) = D_FF(:,:,i)\p_F;                                                     % Displacement's amplitudes with consistent mass matrix

end

% Squeeze of vectors of amplitude in order to simplify graphics

Q0_c = squeeze(q0_c(FDOF,:));

%% FIGURES OF CONSERVATIVE SYSTEM

% Plots Amplitude vs Frecuency ------- Plotted with Y axis as a logarithm
figure(1)
semilogy(f,abs(Q0_c))                   
title("Amplitude Bode Diagram of Conservative System","FontSize",12)
xlabel("Frecuency [Hz]"); ylabel("Amplitude [m]")

% Plots Angular offset vs Frecuency 
figure(2)
plot(f,unwrap(angle(Q0_c)))                  
title("Angular offset Bode Diagram of Conservative System","FontSize",12)
xlabel("Frecuency [Hz]"); ylabel("Phase [rad]")

fprintf('Conservative system finished\n');

%% PROPORTIONAL DAMPING MODEL ([F] = α[M] + β[K])

% Dumping matrix with consistent mass matrix
F_c = a*M_FF + b*K_FF;

% RESOLUTION OF THE NON-CONSERVATIVE DYNAMIC SYSTEM (q0[[K] - Ω^2[M] + i*Ω*[F]] = p0)

D_d_FF = zeros(FDOF,FDOF,F); 
q0_d = zeros(FDOF,F);

for i = 1:F           % This loop makes a sweep in the frecuencies up to the maximum frecuency of interest (Hz)

    D_d_FF(:,:,i) = (K_FF - (2*pi*i)^2*M_FF + 1i*(2*pi*i)*F_c);           % Dynamic stiffness matrix with consistent mass matrix
    q0_d(:,i) = D_d_FF(:,:,i)\p_F;                                        % Displacement's amplitudes with consistent mass matrix

end

% Squeeze of vectors of amplitude in order to simplify graphics

Q0_d = squeeze(q0_d(FDOF,:));

%% FIGURES OF NON-CONSERVATIVE SYSTEM

% Plots Amplitude vs Frecuency ------- Plotted with Y axis as a logarithm
figure(3)
semilogy(f,abs(Q0_d))                
title("Amplitude Bode Diagram of Non-Conservative System","FontSize",12)
xlabel("Frecuency [Hz]"); ylabel("Amplitude [m]");

% Plots Angular offset vs Frecuency 
figure(4)
plot(f,angle(Q0_d))                  
title("Angular offset Bode Diagram of Non-Conservative System","FontSize",12)
xlabel("Frecuency [Hz]"); ylabel("Phase [rad]");

fprintf('Non-Conservative system finished\n');
