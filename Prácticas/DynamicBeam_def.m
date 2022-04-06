% *************************************************************
%           FEM DYNAMIC ANALYSIS OF 1D FREE-FREE BEAM 
% *************************************************************

clc;clear;close all;
set(0,'DefaultFigureVisible','on')

%% MASIC AND GEOMETRIC DATA (all in ISU)

beam.E = 2E9;                                       % Elastic modulus
beam.L = 0.2;                                       % Beam's length
beam.b = 0.02;                                      % Beam's width
beam.t = 0.004;                                     % Beam's thickness
beam.m = 0.03;                                      % Beam's mass
beam.Ixx = beam.b*beam.t^3/12;                      % Beam's area moment of inertia
beam.rho = beam.m/beam.L/beam.b/beam.t;             % Beam's density
beam.rhom = beam.rho*beam.b*beam.t;                 % Beam's linear mass density
beam.etal = [4.73 10.996 17.27876];                 % Free-Free beam's ηl coeficients
beam.x = (0:0.001:beam.L);                          % Beam's partition
c = sqrt(beam.E*beam.Ixx/beam.rho/beam.b/beam.t);   % Constant to be used in continuous model 
a = 1E-6; b = 1E-6;                                 % Propotional damping coefficient

%% INPUT DATA

ne = 100;                       % Number of elements to be used (determined by wavelenght and propagation speed of the wave in the beam)
nn = ne + 1;                    % Number of nodes
dofn = 2;                       % Degrees of freedom per node (only considering flexion)
DOF = dofn*nn;                  % Total dof

p = zeros(DOF,1);
p(DOF/2)= 1;                    % Input force's amplitudes (each value represents deflection and twist of each node of the beam)
F = 2000;                       % Maximum frecuency
f = (1:1:F);                    % Frecuency sweep of the input force

beam.modes = (1:7);             % Shape mode plotted (ascending order: first 2 are rigid solid modes and the other are the ones we want)
beam.resonance_i = [3 5 7];     % Index i is where the resonance frequencies are according to the theory

%% STIFFNESS AND INERTIA BEAM MATRICES

coord_n = zeros(nn,2);                  % Nodal coordinates matrix
coord_n(:,1) = 0:beam.L/ne:beam.L;      % Only X coordinate is different than 0

connect_e = zeros(ne,2);                % Connectivity matrix of elements through nodes
connect_e(:,1) = 1:1:ne;
connect_e(:,2) = 2:1:nn;

K = zeros(DOF);                         % Initilization of the stiffness matrix
M_consist = zeros(DOF);                 % Consistent mass matrix
M_lumped = zeros(DOF);                  % Lumped mass matrix

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

    m_e = beam.rho*beam.b*beam.t*Le;            % Mass of the element          
    Mce = m_e/420*[156 22*Le 54 -13*Le;...      % Consistent mass matrix of each element
          22*Le 4*Le^2 13*Le -3*Le^2;...
          54 13*Le 156 -22*Le;...
          -13*Le -3*Le^2 -22*Le 4*Le^2];     

    M_consist(dofe,dofe) = M_consist(dofe,dofe) + Mce;

 % Second compute the lumped mass matrix

    Mle = m_e/2*[1 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0];       % Lumped mass matrix of each element

    M_lumped(dofe,dofe) = M_lumped(dofe,dofe) + Mle;
end

%% RESOLUTION OF THE CONSERVATIVE DYNAMIC SYSTEM (q0[[K] - Ω^2[M]] = p0)

% Dumped is assumed negligible

D_c = zeros(DOF,DOF,F); D_l = zeros(DOF,DOF,F);
q0_c = zeros(DOF,F); q0_l = zeros(DOF,F);
q_c = zeros(DOF,F); q_l = zeros(DOF,F);

for i = 1:F           % This loop makes a sweep in the frecuencies up to the maximum frecuency of interest (Hz)

    D_c(:,:,i) = (K - (2*pi*i)^2*M_consist);                                              % Dynamic stiffness matrix with consistent mass matrix
    q0_c(:,i) = D_c(:,:,i)\p;                                                             % Displacement's amplitudes with consistent mass matrix

    D_l(:,:,i) = (K - (2*pi*i)^2*M_lumped);                                               % Dynamic stiffness matrix with lumped mass matrix
    q0_l(:,i) = D_l(:,:,i)\p;                                                             % Displacement's amplitudes with lumped mass matrix

end

% Natural frecuencies (to be compared with graphics)

w1 = 0;                                        % Rigid-body motion
w2 = beam.etal(1,1)^2*c/beam.L^2/2/pi;
w3 = beam.etal(1,2)^2*c/beam.L^2/2/pi;
w4 = beam.etal(1,3)^2*c/beam.L^2/2/pi;

% Squeeze of vectors of amplitude in order to simplify graphics

Q0_c = squeeze(q0_c(DOF/2,:));
Q0_l = squeeze(q0_l(DOF/2,:));

%% FIGURES OF CONSERVATIVE SYSTEM

% Plots Amplitude vs Frecuency ------- Plotted with Y axis as a logarithm
fig1 = figure;
semilogy(f,abs(Q0_c))                
hold on
semilogy(f,abs(Q0_l))      
title("Amplitude Bode Diagram of Conservative System","FontSize",12)
legend("Consistent mass matrix", "Lumped mass matrix")
xlabel("Frecuency [Hz]"); ylabel("Amplitude [m]")
savefig(fig1,'fig1.fig')

% Find the resonance and antiresonance frequencies
[peaks.rcons_c,locs_rc] = findpeaks(abs(Q0_c));
[peaks.rcons_l,locs_rl] = findpeaks(abs(Q0_l));
text(locs_rl+.02,peaks.rcons_l,num2str((1:numel(peaks.rcons_l))'))

[peaks.acons_c,locs_ac] = findpeaks(abs(1./Q0_c));
[peaks.acons_l,locs_al] = findpeaks(abs(1./Q0_l));

% Plots Angular offset vs Frecuency 
fig2 = figure;
plot(f,unwrap(angle(Q0_c)))                  
hold on
plot(f,unwrap(angle(Q0_l)))
title("Angular offset Bode Diagram of Conservative System","FontSize",12)
legend("Consistent mass matrix", "Lumped mass matrix")
xlabel("Frecuency [Hz]"); ylabel("Phase [rad]")
savefig(fig2,'fig2.fig')

fprintf('Conservative system finished\n');

%% PROPORTIONAL DAMPING MODEL ([F] = α[M] + β[K])

% Dumping matrix with consistent mass matrix
F_c = a*M_consist + b*K;
% Dumping matrix with lumped mass matrix
F_l = a*M_lumped + b*K;

%% RESOLUTION OF THE NON-CONSERVATIVE DYNAMIC SYSTEM (q0[[K] - Ω^2[M] + i*Ω*[F]] = p0)

D_dc = zeros(DOF,DOF,F); D_dl = zeros(DOF,DOF,F);
q0_dc = zeros(DOF,F); q0_dl = zeros(DOF,F);
q_dc = zeros(DOF,F); q_dl = zeros(DOF,F);

for i = 1:F           % This loop makes a sweep in the frecuencies up to the maximum frecuency of interest (Hz)

    D_dc(:,:,i) = (K - (2*pi*i)^2*M_consist + 1i*(2*pi*i)*F_c);           % Dynamic stiffness matrix with consistent mass matrix
    q0_dc(:,i) = D_dc(:,:,i)\p;                                           % Displacement's amplitudes with consistent mass matrix

    D_dl(:,:,i) = (K - (2*pi*i)^2*M_lumped + 1i*(2*pi*i)*F_c);            % Dynamic stiffness matrix with lumped mass matrix
    q0_dl(:,i) = D_dl(:,:,i)\p;                                           % Displacement's amplitudes with lumped mass matrix

end

% Squeeze of vectors of amplitude in order to simplify graphics

Q0_dc = squeeze(q0_dc(DOF/2,:));
Q0_dl = squeeze(q0_dl(DOF/2,:));

%% FIGURES OF NON-CONSERVATIVE SYSTEM

% Plots Amplitude vs Frecuency ------- Plotted with Y axis as a logarithm
fig3 = figure;
semilogy(f,abs(Q0_dc))                
title("Amplitude Bode Diagram of Non-Conservative System","FontSize",12)
xlabel("Frecuency [Hz]"); ylabel("Amplitude [m]");
hold on
semilogy(f,abs(Q0_dl))                       
legend("Consistent mass matrix", "Lumped mass matrix")
savefig(fig3,'fig3.fig')

% Find the resonance and antiresonance frequencies
[peaks.rnocons_c,locs_nocons_rc] = findpeaks(abs(Q0_dc));
[peaks.rnocons_l,locs_nocons_rl] = findpeaks(abs(Q0_dl));
text(locs_nocons_rl+.02,peaks.rnocons_l,num2str((1:numel(peaks.rnocons_l))'))

[peaks.acons_l,locs_nocons_ac] = findpeaks(abs(1./Q0_dc));
[peaks.acons_l,locs_nocons_al] = findpeaks(abs(1./Q0_dl));

%disp("Antiresonance frequencies obteained with the Bode diagram")


% Plots Angular offset vs Frecuency 
fig4 = figure;
plot(f,angle(Q0_dc))                  
title("Angular offset Bode Diagram of Non-Conservative System","FontSize",12)
xlabel("Frecuency [Hz]"); ylabel("Phase [rad]");
hold on
plot(f,angle(Q0_dl))
legend("Consistent mass matrix", "Lumped mass matrix")
savefig(fig4,'fig4.fig')

fprintf('Non-Conservative system finished\n');

%% MODAL SHAPES OF DISCREET MODEL (Natural frequencies are the eigenvalues of the problem and modal shapes are the eigenvectors)

% With consistent mass matrix
[V_consist,D_consist] = eig(K,M_consist);
[W_c,order_c] = sort(sum(sqrt(D_consist)./2./pi));                          % Vector of natural frequencies ordered in ascending order
W_c = diag(W_c);                                                            % Natural frequencies in a diagonal matrix
V_consist = V_consist(:,order_c);                                           % Rewrite modal shapes so they match the adcending order of frequencies
% With lumped mass matrix
[V_lumped,D_lumped] = eig(K,M_lumped);
[W_l,order_l] = sort(sum(sqrt(D_lumped)./2./pi,1));                         % Vector of natural frequencies ordered in ascending order
W_l = diag(W_l);                                                            % Natural frequencies in a diagonal matrix
V_lumped = V_lumped(:,order_l);

disp("Frequencies obtained by solving the eigenvalue problem")

vw_eig_c = [W_c(beam.resonance_i(:,1),beam.resonance_i(:,1)),...
            W_c(beam.resonance_i(:,2),beam.resonance_i(:,2)),...
            W_c(beam.resonance_i(:,3),beam.resonance_i(:,3))]';

vw_eig_l = [W_l(beam.resonance_i(:,1),beam.resonance_i(:,1)),...
            W_l(beam.resonance_i(:,2),beam.resonance_i(:,2)),...
            W_l(beam.resonance_i(:,3),beam.resonance_i(:,3))]';             % Starts in 3rd because 1 and 2 are rigid solid modes 

disp(num2str([vw_eig_c, vw_eig_l]))

disp("Antiresonance frequencies obtained through the graphics")
disp(num2str(locs_ac'))

x_axis = linspace(0,beam.L,DOF/2);

fig5 = figure;
plot(x_axis,V_consist(1:2:DOF,beam.modes(:)))
set(gca,'YTick',[])
title("Mode shapes of discreet model with consistent mass matrix")
xlabel("X-coordinate [m]"); ylabel("Deformation")
legend("First Mode", "Second Mode", "Third Mode")
savefig(fig5,'fig5.fig')

fig6 = figure;
plot(x_axis,V_lumped(1:2:DOF,beam.modes(3:5)))
set(gca,'YTick',[])
title("Mode shapes of discreet model with lumped mass matrix")
xlabel("X-coordinate [m]"); ylabel("Deformation")
legend("First Mode", "Second Mode", "Third Mode")
savefig(fig6,'fig6.fig')

% Compute kj and mj with the discreet model in the modal space

mj_c = sum(V_consist(:,beam.modes)'*M_consist*V_consist(:,beam.modes));
kj_c = sum(V_consist(:,beam.modes)'*K*V_consist(:,beam.modes));

mj_l = sum(V_lumped(:,beam.modes)'*M_lumped*V_lumped(:,beam.modes));
kj_l = sum(V_lumped(:,beam.modes)'*K*V_lumped(:,beam.modes));

% Compute ωj with kj and mj

w_c = sqrt(kj_c./mj_c)/2/pi;
w_l = sqrt(kj_l./mj_l)/2/pi;

disp("Frequencies with mj and kj in modal space")

vw_modsp_c = [w_c(:,beam.resonance_i(:,1)),...              % The frequencies asociated are in the positions 3, 5 and 7
              w_c(:,beam.resonance_i(:,2)),...
              w_c(:,beam.resonance_i(:,3))]';

vw_modsp_l = [w_l(:,beam.resonance_i(:,1)),...
              w_l(:,beam.resonance_i(:,2)),...
              w_l(:,beam.resonance_i(:,3))]';

disp(num2str([vw_modsp_c, vw_modsp_l]))

fprintf('Modal shapes of the discreet model finished\n');

%% CONTINUOUS MODEL (Theory of vibration vol II (4.3), Shabana)

% Mode shapes (Ф_j(x))

D = zeros(size(beam.etal));
phi = zeros(size(beam.etal,2),size(beam.x,2));

for i = 1:size(beam.etal,2)
    D(:,i) = - ((cosh(beam.etal(i)) - cos(beam.etal(i)))/(sinh(beam.etal(i)) + sin(beam.etal(i))));
    phi(i,:) = -(sinh(beam.etal(i)/beam.L*beam.x) + sin(beam.etal(i)/beam.L*beam.x) + D(i)*(cosh(beam.etal(i)/beam.L*beam.x) + cos(beam.etal(i)/beam.L*beam.x)));
end

fig7 = figure;
plot(beam.x,phi)
set(gca,'YTick',[])
title("First three mode shapes of a beam with free free ends")
legend("First mode", "Second mode", "Third mode")
xlabel("X-coordinate [m]"); ylabel("Deformation")
savefig(fig7,'fig7.fig')

% Time response (q(t))

mj = zeros(1,size(beam.etal,2));
kj = zeros(1,size(beam.etal,2));
q0 = zeros(size(beam.etal,2),size(f,2));

for i = 1:size(beam.etal,2)
    mj(i) = beam.rho*beam.b*beam.t*trapz(beam.x,phi(i,:).^2,2);                             % Equivalent mass
    kj(i) = beam.E*beam.Ixx*trapz(beam.x(1:end-2),diff(phi(i,:),2).^2,2)*1E12;              % Equivalent stiffness                 
    q0(i,:) = p(DOF/2)*phi(i,101)./(-mj(i)*(2*pi*f(:)).^2 + kj(i));                         % Modal coordinates
end

% Beam's transverse vibration (v(x,t) = (ΣФ(x)*q0(Ω))*exp(iΩ*t))

v0 = phi'*q0;                   % Amplitude of the vibration

fig8 = figure;
semilogy(f,abs(v0(101,:)))
title("Response of a continuous free-free beam")
xlabel("Frecuency [Hz]"); ylabel("Transverse vibration [m]")
savefig(fig8,'fig8.fig')

% Find the resonance and antiresonance frequencies
[peaks.rcont,locsr_cont] = findpeaks(abs(v0(101,:)));
[peaks.acont,locsa_cont] = findpeaks(abs(1./v0(101,:)));

disp('Frequencies of the discreete and continous model:')
vw_cont = [w2,w3,w4]';
disp(num2str([vw_cont, (sqrt(kj./mj)/2/pi)']))

fprintf('Continuous model finished\n');

%% SIMULATED TEST

promt.ask = 'Want to realize a simulated test? (1/0) ';
answer = input(promt.ask);

if answer == 1
    % Input of frequencies obtained in the dynamic test
    promt.af = 'First three resonance frequencies (Hz): ';
    beam.rf = input(promt.af);
    promt.rf = 'First three antiresonance frequencies (Hz): ';
    beam.af = input(promt.rf);

    % Compute Young's Modulus with 3 different methods

    beam.E_iso = E_ISO(beam.L,beam.rhom,beam.Ixx,beam.af);      % ISO-16940
    beam.E_aff = E_AFF(beam.L,beam.rhom,beam.Ixx,beam.af);      % Antiresonances of a free-free beam
    beam.E_rff = E_RFF(beam.L,beam.rhom,beam.Ixx,beam.rf);      % Resonances of a free-free beam
else
    fprintf('No simulated test conducted\n')
end

%% EXPORT PLOTS AS .TXT DATA FILES

h = openfig('fig8.fig');
h = findobj(gca,'Type','line');
x = get (h,'Xdata');
y = get(h,'Ydata');

A = [];
A(:,1) = x;
A(:,2) = y;

writematrix(A,'F90001ZZa.txt','Delimiter',',')

%% FUNCTIONS

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
