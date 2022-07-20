% ***********************************************************************
%           ANALYSIS OF A CANTILEVER BEAM AS A TUNNED MASS DUMPER 
% ***********************************************************************

 clc;clear;close all;

%% MASIC AND GEOMETRIC DATA (all in ISU)

beam.E = 3E9;                                       % Elastic modulus
beam.L = 0.05;                                       % Beam's length
beam.b = 0.01;                                      % Beam's width
beam.t = 0.001;                                     % Beam's thickness
beam.tl = 0.012;                                     % Thickness of the last element of the beam      
beam.Ixx = beam.b*beam.t^3/12;                      % Beam's area moment of inertia
beam.rho = 1.1678e+03;             % Beam's density
beam.m = beam.rho*beam.L*beam.b*beam.t;
beam.rhom = beam.rho*beam.b*beam.t;                 % Beam's linear mass density
beam.x = (0:0.001:beam.L);                          % Beam's partition 
a = 1E-6; b = 1E-6;                                 % Proportional damping coefficient

%% INPUT DATA

ne = 10;                        % Number of elements to be used (determined by wavelenght and propagation speed of the wave in the beam)
nn = ne + 1;                    % Number of nodes
dofn = 2;                       % Degrees of freedom per node (only considering flexion)
DOF = dofn*nn;                  % Total dof 

tmd_elements = [7 8 9 10];

p = zeros(DOF,1);
p(17)= - 9.81*0.00233356;                    % Input force's amplitudes (each value represents deflection and twist of each node of the beam)
F = 2000;                       % Maximum frecuency
f = (1:1:F);                    % Frecuency sweep of the input force

restr_dof = [1 2];              % Restricted DOF in the fixed end
free_dof = (3:DOF);             % Free DOF                    

RDOF = length(restr_dof);       % Number of restricted DOF
FDOF = DOF - RDOF;              % Number of free DOF

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
            
    if ismember(e,tmd_elements)                                                  % Changes the thickness of the beam in order to simulate a TMD
        beam.t = beam.tl;
        beam.Ixx = beam.b*beam.t^3/12;
    else
        beam.t = beam.t;
    end

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


%% STATIC SYSTEM

u = zeros(dofn*nn - 1,1);

u_F = K_FF\p_F;
u(free_dof,1) = u_F;

w = u(1:2:DOF);

figure()
plot(w)
%% RESOLUTION OF THE CONSERVATIVE DYNAMIC SYSTEM (q0[[K] - Ω^2[M]] = p0)

% Dumped is assumed negligible

D_FF = zeros(FDOF,FDOF,F);
q0_cF = zeros(FDOF,F);
q0_c = zeros(DOF,F);

for i = 1:F           % This loop makes a sweep in the frecuencies up to the maximum frecuency of interest (Hz)

    D_FF(:,:,i) = (K_FF - (2*pi*i)^2*M_FF);                                           % Dynamic stiffness matrix with consistent mass matrix
    q0_cF(:,i) = D_FF(:,:,i)\p_F;                                                     % Displacement's amplitudes with consistent mass matrix
    q0_c(free_dof,i) = q0_cF(:,i);

end

% Squeeze of vectors of amplitude in order to simplify graphics

Q0_c = squeeze(q0_c(FDOF-1,:));

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
q0_dF = zeros(FDOF,F);
q0_d = zeros(DOF,F);

for i = 1:F           % This loop makes a sweep in the frecuencies up to the maximum frecuency of interest (Hz)

    D_d_FF(:,:,i) = (K_FF - (2*pi*i)^2*M_FF + 1i*(2*pi*i)*F_c);           % Dynamic stiffness matrix with consistent mass matrix
    q0_dF(:,i) = D_d_FF(:,:,i)\p_F;                                        % Displacement's amplitudes with consistent mass matrix
    q0_d(free_dof,i) = q0_dF(:,i); 

end

% Squeeze of vectors of amplitude in order to simplify graphics

Q0_d = squeeze(q0_d(DOF,:));

%% NATURAL FREQUENCIES EXPECTED (Bleving Beams)

lambda(1) = 1.875104; lambda(2) = 4.694091; lambda(3) = 7.854757;
w = zeros(3,1);

for i = 1:3
    w(i,1) = lambda(i)^2*sqrt(beam.E*beam.Ixx/beam.rho/beam.b/beam.t)/(2*pi*beam.L^2);          % Natural frequencies that should appear
end


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
