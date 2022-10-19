% *************************************************************
%            FEM DYNAMIC ANALYSIS OF 1D FREE-FREE BEAM 
% *************************************************************

% This file works as a compendium where I have been coding
% different methods and alternatives to do the same things, 
% for the final and optimize code go to aDynamicBeam_def.m

clc;clear;close all;

tic
%% MASIC AND GEOMETRIC DATA (all in ISU)

beam.E = 2E9;                                           % Elastic modulus
beam.L = 0.2;                                           % Beam's length
beam.b = 0.02;                                          % Beam's width
beam.t = 0.004;                                         % Beam's thickness
beam.Ixx = beam.b*beam.t^3/12;                          % Beam's area moment of inertia
beam.m = 0.03;                                          % Beam's mass
beam.rho = beam.m/beam.L/beam.b/beam.t;                 % Beam's density
beam.rhom = beam.rho*beam.b*beam.t;                     % Beam's linear mass density

beam.etal = [4.73 10.996 17.27876];                     % Free-Free beam's ηl coeficients
beam.sigma = [0.982502215 1.000777312 0.99996645];      % Sigma coeficients
beam.x = (0:0.001:beam.L);                              % Beam's partition

promt.af = 'First three resonance frequencies (Hz): ';
beam.rf = input(promt.af);
promt.rf = 'First three antiresonance frequencies (Hz): ';
beam.af = input(promt.rf);

beam.E_iso = E_ISO(beam.L,beam.rhom,beam.Ixx,beam.af);
beam.E_aff = E_AFF(beam.L,beam.rhom,beam.Ixx,beam.af);
beam.E_rff = E_RFF(beam.L,beam.rhom,beam.Ixx,beam.rf);

%% INPUT DATA

ne = 20;                        % Number of elements to be used (determined by wavelenght and propagation speed of the wave in the beam)
nn = ne + 1;                    % Number of nodes
dofn = 2;                       % Degrees of freedom per node (only considering flexion)
DOF = dofn*nn;                  % Total dof

p = zeros(DOF,1);               % Input force's amplitudes (each value represents deflection and twist of each node of the beam)
p(3)= 1;                        
F = 2000;                       % Maximum frecuency
f = (1:1:F);                    % Frecuency sweep of the input force

mode = 3;                       % DOF plotted
beam.modes = (1:5);             % Shape mode plotted (ascending order: first 2 are rigid solid and the other 3 are the ones we want)

c = sqrt(beam.E*beam.Ixx/beam.rho/beam.t/beam.b);       % Constant to be used in continuous model

a = 1E-4; b = 1E-4;                  % Propotional damping coefficient

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
chi_ip.ff = [-1/3*sqrt(5+2*sqrt(10/7))  -1/3*sqrt(5-2*sqrt(10/7)) 0 1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7))];   
w_ip.ff = [(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900]; 

%% STIFFNESS BEAM MATRIX

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
        chi = chi_ip.ff(i);
        w = w_ip.ff(i);
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

    Mle = beam.rho*beam.b*beam.t*Le/2*[1 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0];

    M_lumped(dofe,dofe) = M_lumped(dofe,dofe) + Mle;
end

%% RESOLUTION OF THE CONSERVATIVE DYNAMIC SYSTEM (q0[[K] - Ω^2[M]] = p0)
% Dumped is assumed negligible

D_c = zeros(DOF,DOF,F); D_l = zeros(DOF,DOF,F);
q0_c = zeros(DOF,F); q0_l = zeros(DOF,F);
q_c = zeros(DOF,F); q_l = zeros(DOF,F);

for i = 1:F           % This loop makes a sweep in the frecuencies up to the maximum frecuency of interest (Hz)

    D_c(:,:,i) = (K - (2*pi*i)^2*M_consist_2);         % Dynamic stiffness matrix with consistent mass matrix
    q0_c(:,i) = D_c(:,:,i)\p;                          % Displacement's amplitudes with consistent mass matrix

    D_l(:,:,i) = (K - (2*pi*i)^2*M_lumped);            % Dynamic stiffness matrix with lumped mass matrix
    q0_l(:,i) = D_l(:,:,i)\p;                          % Displacement's amplitudes with lumped mass matrix

end

% Natural frecuencies (to be compared with graphics)

w1 = 0;                             % Rigid-body motion
w2 = 4.730^2*c/beam.L^2/2/pi;
w3 = 7.853^2*c/beam.L^2/2/pi;
w4 = 10.996^2*c/beam.L^2/2/pi;

% Squeeze of vectors of amplitude in order to simplify graphics

Q0_c = squeeze(q0_c(mode,:));
Q0_l = squeeze(q0_l(mode,:));

%% FIGURES OF CONSERVATIVE SYSTEM
% Plots Amplitude vs Frecuency ------- Plotted with Y axis as a logarithm
close all
figure(1)
semilogy(f,abs(Q0_c))                
title("Amplitude Bode Diagram of Conservative System","FontSize",12)
xlabel("Frecuency [Hz]"); ylabel("Amplitude [m]");
hold on
semilogy(f,abs(Q0_l))                       
legend("Consistent mass matrix", "Lumped mass matrix")
 
% Plots Angular offset vs Frecuency 
figure(2)
plot(f,unwrap(angle(Q0_c)))                  
title("Angular offset Bode Diagram of Conservative System","FontSize",12)
xlabel("Frecuency [Hz]"); ylabel("Phase [rad]");
hold on
plot(f,unwrap(angle(Q0_l)))
legend("Consistent mass matrix", "Lumped mass matrix")

%% PROPORTIONAL DAMPING MODEL ([F] = α[M] + β[K])

% Dumping matrix with consistent mass matrix
F_c = a*M_consist_2 + b*K;
% Dumping matrix with lumped mass matrix
F_l = a*M_lumped + b*K;

%% RESOLUTION OF THE NON-CONSERVATIVE DYNAMIC SYSTEM (q0[[K] - Ω^2[M] + i*Ω*[F]] = p0)

D_dc = zeros(DOF,DOF,F); D_dl = zeros(DOF,DOF,F);
q0_dc = zeros(DOF,F); q0_dl = zeros(DOF,F);
q_dc = zeros(DOF,F); q_dl = zeros(DOF,F);

for i = 1:F           % This loop makes a sweep in the frecuencies up to the maximum frecuency of interest (Hz)

    D_dc(:,:,i) = (K - (2*pi*i)^2*M_consist_2 + 1i*(2*pi*i)*F_c);         % Dynamic stiffness matrix with consistent mass matrix
    q0_dc(:,i) = D_dc(:,:,i)\p;                                           % Displacement's amplitudes with consistent mass matrix

    D_dl(:,:,i) = (K - (2*pi*i)^2*M_lumped + 1i*(2*pi*i)*F_c);            % Dynamic stiffness matrix with lumped mass matrix
    q0_dl(:,i) = D_dl(:,:,i)\p;                                           % Displacement's amplitudes with lumped mass matrix

end

% Squeeze of vectors of amplitude in order to simplify graphics

Q0_dc = squeeze(q0_dc(mode,:));
Q0_dl = squeeze(q0_dl(mode,:));

%% FIGURES OF NON-CONSERVATIVE SYSTEM
% Plots Amplitude vs Frecuency ------- Plotted with Y axis as a logarithm
figure(3)
semilogy(f,abs(Q0_dc))                
title("Amplitude Bode Diagram of Non-Conservative System","FontSize",12)
xlabel("Frecuency [Hz]"); ylabel("Amplitude [m]");
hold on
semilogy(f,abs(Q0_dl))                       
legend("Consistent mass matrix", "Lumped mass matrix")
 
% Plots Angular offset vs Frecuency 
figure(4)
plot(f,angle(Q0_dc))                  
title("Angular offset Bode Diagram of Non-Conservative System","FontSize",12)
xlabel("Frecuency [Hz]"); ylabel("Phase [rad]");
hold on
plot(f,angle(Q0_dl))
legend("Consistent mass matrix", "Lumped mass matrix")

%% MODAL SHAPES OF DISCREET MODEL (Natural frequencies are the eigenvalues of the problem and modal shapes are the eigenvectors)

% With consistent mass matrix
[V_consist,D_consist] = eig(K,M_consist_2);
[W_c,order_c] = sort(sum(sqrt(D_consist)./2./pi));                          % Vector of natural frequencies ordered in ascending order
W_c = diag(W_c);                                                            % Natural frecuencies in a diagonal matrix
V_consist = V_consist(:,order_c);                                           % Rewrite modal shapes so they match the adcending order of frequencies
% With lumped mass matrix
[V_lumped,D_lumped] = eig(K,M_lumped);
[W_l,order_l] = sort(sum(sqrt(D_lumped)./2./pi,1));                         % Vector of natural frequencies ordered in ascending order
W_l = diag(W_l);                                                            % Natural frecuencies in a diagonal matrix
V_lumped = V_lumped(:,order_l);

disp("Frequencies obtained by solving the eigenvalue problem")
vw_eig_c = [W_c(3,3), W_c(4,4), W_c(5,5)]'; vw_eig_l = [W_l(3,3), W_l(4,4), W_l(5,5)]';   % Starts in 3rd because 1 and 2 are rigid solid modes 
disp(num2str([vw_eig_c, vw_eig_l]))

x = linspace(0,beam.L,DOF/2);

figure(5)
plot(x,V_consist(1:2:DOF,beam.modes(3:5)))
set(gca,'YTick',[])
title("Mode shapes of discreet model with consistent mass matrix")
xlabel("X-coordinate [m]"); ylabel("Deformation")
legend("First Mode", "Second Mode", "Third Mode")

figure(6)
plot(x,V_lumped(1:2:DOF,beam.modes(3:5)))
set(gca,'YTick',[])
title("Mode shapes of discreet model with lumped mass matrix")
xlabel("X-coordinate [m]"); ylabel("Deformation")
legend("First Mode", "Second Mode", "Third Mode")

% Compute kj and mj with the discreet model in the modal space

mj_c = sum(V_consist(:,beam.modes)'*M_consist*V_consist(:,beam.modes));
kj_c = sum(V_consist(:,beam.modes)'*K*V_consist(:,beam.modes));

mj_l = sum(V_lumped(:,beam.modes)'*M_lumped*V_lumped(:,beam.modes));
kj_l = sum(V_lumped(:,beam.modes)'*K*V_lumped(:,beam.modes));

% Compute ωj with kj and mj

w_c = sqrt(kj_c./mj_c)/2/pi;
w_l = sqrt(kj_l./mj_l)/2/pi;

disp("Frequencies with mj and kj in modal space")
vw_modsp_c = [w_c(1,3),w_c(1,4), w_c(1,5)]'; vw_modsp_l = [w_l(1,3), w_l(1,4), w_l(1,5)]';
disp(num2str([vw_modsp_c, vw_modsp_l]))

%% CONTINUOUS MODEL (Theory of vibration vol II (4.3), Shabana)

% Mode shapes (Ф_j)

    % Done with Shabana equation

D = zeros(size(beam.etal));
phi = zeros(size(beam.etal,2),size(beam.x,2));

for i = 1:size(beam.etal,2)
    D(:,i) = - ((cosh(beam.etal(i)) - cos(beam.etal(i)))/(sinh(beam.etal(i)) + sin(beam.etal(i))));
    phi(i,:) = -(sinh(beam.etal(i)/beam.L*beam.x) + sin(beam.etal(i)/beam.L*beam.x) + D(i)*(cosh(beam.etal(i)/beam.L*beam.x) + cos(beam.etal(i)/beam.L*beam.x)));
end

figure(7)
plot(beam.x,phi)
title("First three mode shapes of a beam with free free ends")
legend("First mode", "Second mode", "Third mode")
xlabel("X-coordinate [m]"); ylabel("Deformation")

    % Done with Blevins sigma coeficient

for i = 1:size(beam.etal,2)
    phi(i,:) = -beam.sigma(i)*(sinh(beam.etal(i)/beam.L*beam.x) + sin(beam.etal(i)/beam.L*beam.x)) + (cosh(beam.etal(i)/beam.L*beam.x) + cos(beam.etal(i)/beam.L*beam.x));
end

figure(8)
plot(beam.x,phi)
title("First three mode shapes of a beam with free free ends")
legend("First mode", "Second mode", "Third mode")
xlabel("X-coordinate [m]"); ylabel("Deformation")

% Figures 5 and 6 should be the same

% Time response (q(t))

mj = zeros(1,size(beam.etal,2));
kj = zeros(1,size(beam.etal,2));
q0 = zeros(size(beam.etal,2),size(f,2));

mj_ = zeros(1,size(beam.etal,2));
kj_ = zeros(1,size(beam.etal,2));
q0_ = zeros(size(beam.etal,2),size(f,2));
 
for i = 1:size(beam.etal,2)             % Done with trapz function
    mj_(i) = beam.rho*beam.b*beam.t*trapz(beam.x,phi(i,:).^2,2);                     % Equivalent mass
    kj_(i) = beam.E*beam.Ixx*trapz(beam.x(1:end-2),diff(phi(i,:),2).^2,2)*1E12;           % Equivalent stiffness
    q0_(i,:) = p(3)*phi(i,101)./(-mj_(i)*(2*pi*f(:)).^2 + kj_(i));                     % Modal coordinates
end

for i = 1:size(beam.etal,2)             % Done with Simpson's integration rule
    mj(i) = beam.rho*beam.b*beam.t*simps(beam.x,phi(i,:).^2,2);                      % Equivalent mass
    kj(i) = beam.E*beam.Ixx*simps(beam.x(1:end-2),(diff(phi(i,:),2)).^2,2)*1E12;          % Equivalent stiffness
    q0(i,:) = p(3)*phi(i,101)./(-mj(i)*(2*pi*f(:)).^2 + kj(i));                      % Modal coordinates
end

% Beam's transverse vibration (v(x,t) = (ΣФ(x)*q0(Ω))*exp(iΩ*t))

v0 = phi'*q0;                   % Amplitude of the vibration
v0_ = phi'*q0_;

figure(9)
semilogy(f,abs(abs(v0(:,:))))
title("Response of a continuous free-free beam using Simpson's rule")
xlabel("Frecuency [Hz]"); ylabel("Transverse vibration [m]")

figure(10)
semilogy(f,abs(abs(v0(:,:))))
title("Response of a continuous free-free beam unsing 'trapz' function")
xlabel("Frecuency [Hz]"); ylabel("Transverse vibration [m]")

disp('Frequencies of the discreete and continous model (Simpsons rule and trapz)')
vw = [w2,w3,w4]';
disp(num2str([vw, (sqrt(kj./mj)/2/pi())', (sqrt(kj_./mj_)/2/pi)']))

toc
%% FUNCTIONS

function NH = Nh(chi,Je)          % Hermitic form functions
    NH = 1/4*[(2-3*chi+chi^3) Je*(1-chi-chi^2+chi^3) (2+3*chi-chi^3) Je*(-1-chi+chi^2+chi^3)];
end

function sd_NH = sd_Nh(chi,Je)    % Second derivative of hermitic form functions
    sd_NH = [(chi*3/2) Je*(chi*3-1)/2 -(chi*3/2) Je*(chi*3+1)/2];
end

function E_I = E_ISO(L,rhom,Ixx,f)                           % Young's modulus with ISO-6940 method
    lambda = [1.87510 4.69410 7.85476];
    E_I = rhom/Ixx*(2*pi*(L/2)^2*f./lambda.^2).^2*1E-9;
end

function E_A = E_AFF(L,rhom,Ixx,f)                           % Young's modulus with antiresonance frequencies of a free-free beam method
    lambda = [3.75038 9.39740 15.73438];
    E_A = rhom/Ixx*(2*pi*L^2*f./lambda.^2).^2*1E-9;
end

function E_R = E_RFF(L,rhom,Ixx,f)                           % Young's modulus with resonance frequencies of a free-free beam method
    lambda = [4.73004 10.99561 17.27876];
    E_R = rhom/Ixx*(2*pi*L^2*f./lambda.^2).^2*1E-9;
end

