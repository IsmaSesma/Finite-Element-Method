% *************************************************************
%            FEM DYNAMIC ANALYSIS OF 1D FREE-FREE BEAM 
% *************************************************************

clc;clear;close all;

tic
%% MASIC AND GEOMETRIC DATA (all in ISU)

beam.E = 2E9;                                       % Elastic modulus
beam.L = 0.2;                                       % Beam's length
beam.b = 0.02;                                      % Beam's width
beam.t = 0.004;                                     % Beam's thickness
beam.Ixx = beam.b*beam.t^3/12;                      % Beam's area moment of inertia
beam.m = 0.03;                                      % Beam's ,ass
beam.rho = beam.m/beam.L/beam.b/beam.t;             % Beam's density
beam.eta = [4.73 7.853 10.996];                     % Free-Free beam's η coeficients
beam.x = [0:0.001:beam.L];                          % Beam's partition
c = sqrt(beam.E*beam.Ixx/beam.rho/beam.b/beam.t);   % Constant to be used in continuous model   

%% INPUT DATA

ne = 2;                         % Number of elements to be used (determined by wavelenght and propagation speed of the wave in the beam)
nn = ne + 1;                    % Number of nodes
dofn = 2;                       % Degrees of freedom per node (only considering flexion)
DOF = dofn*nn;                  % Total dof

p = zeros(DOF,1);
p(3)= 1;                        % Input force's amplitudes (each value represents deflection and twist of each node of the beam)
F = 800;                        % Maximum frecuency
f = (1:1:F);                    % Frecuency sweep of the input force
mode = 3;                       % DOF plotted

%% STIFFNESS AND INERTIA BEAM MATRICES


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

    m = beam.rho*beam.b*beam.t*Le;            % Mass of the element          
    Mce = m/420*[156 22*Le 54 -13*Le; 22*Le 4*Le^2 13*Le -3*Le^2; 54 13*Le 156 -22*Le; -13*Le -3*Le^2 -22*Le 4*Le^2];

    M_consist(dofe,dofe) = M_consist(dofe,dofe) + Mce;

 % Second compute the lumped mass matrix

    Mle = m/2*[1 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0];

    M_lumped(dofe,dofe) = M_lumped(dofe,dofe) + Mle;
end

%% RESOLUTION OF THE DYNAMIC SYSTEM (q0[[K] - w^2[M]] = p0)
% Dumped is assumed negligible

D_c = zeros(DOF,DOF,F); D_l = zeros(DOF,DOF,F);
q0_c = zeros(DOF,F); q0_l = zeros(DOF,F);
q_c = zeros(DOF,F); q_l = zeros(DOF,F);

for i = 1:F           % This loop makes a sweep in the frecuencies up to the maximum frecuency of interest (Hz)

    D_c(:,:,i) = (K - (2*pi*i)^2*M_consist);           % Dynamic stiffness matrix with consistent mass matrix
    q0_c(:,i) = D_c(:,:,i)\p;                          % Displacement's amplitudes with consistent mass matrix

    D_l(:,:,i) = (K - (2*pi*i)^2*M_lumped);            % Dynamic stiffness matrix with lumped mass matrix
    q0_l(:,i) = D_l(:,:,i)\p;                          % Displacement's amplitudes with lumped mass matrix

    q_c(:,i) = q0_c(:,i)*exp(2*pi*1i*i);               % Complex displacement with consistent mass matrix
    q_l(:,i) = q0_l(:,i)*exp(2*pi*1i*i);               % Complex displacement with lumped mass matrix

end

% Natural frecuencies (to be compared with graphics)

w1 = 0;                     % Rigid-body motion
w2 = 4.730^2*c/beam.L^2/2/pi;
w3 = 7.853^2*c/beam.L^2/2/pi;
w4 = 10.996^2*c/beam.L^2/2/pi;

% Squeeze of vectors of amplitude in order to simplify graphics

Q0_c = squeeze(q0_c(mode,:));
Q0_l = squeeze(q0_l(mode,:));
Q_c = squeeze(q_c(mode,:));
Q_l = squeeze(q_l(mode,:));

%% FIGURES

% Plots Amplitude vs Frecuency ------- Plotted with Y axis as a logarithm
close all
figure(1)
semilogy(f,abs(Q0_c))                
title("Amplitude Bode Diagram","FontSize",12)
hold on
semilogy(f,abs(Q0_l))                       
legend("Consistent mass matrix", "Lumped mass matrix")
xlabel("Frecuency [Hz]"); ylabel("Amplitude [m]")

% Plots Angular offset vs Frecuency 
figure(2)
plot(f,unwrap(angle(Q_c)))                  
title("Angular offset Bode Diagram","FontSize",12)
hold on
plot(f,unwrap(angle(Q_l)))
legend("Consistent mass matrix", "Lumped mass matrix")
xlabel("Frecuency [Hz]"); ylabel("Phase [rad]")

%% CONTINUOUS MODEL (Theory of vibration vol II (4.3), Shabana)

% Mode shapes (Ф_j)

D = zeros(size(beam.eta));
phi = zeros(size(beam.eta,2),size(beam.x,2));

for i = 1:size(beam.eta,2)
    D(:,i) = - (cosh(beam.eta(i)*beam.L) - cos(beam.eta(i)*beam.L)/(sinh(beam.eta(i)*beam.L) + sin(beam.eta(i)*beam.L)));
    phi(i,:) = sinh(beam.eta(i)*beam.x) + sin(beam.eta(i)*beam.x) + D(i)*(cosh(beam.eta(i)*beam.x) + cos(beam.eta(i)*beam.x));
end

figure(3)
plot(beam.x,phi)
title("First four mode shapes of a beam with free free ends")
legend("First mode", "Second mode", "Third mode")
xlabel("X-coordinate [m]"); ylabel("Deformation")

toc




