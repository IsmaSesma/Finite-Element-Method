%% CODE TO PRODUCE GRAPHS FOR THE TFG

clc, close, clear

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',2)

%% ONE-DEGREE OF FREEDOM 

% Variables
K = 1; K_a = 1; F = 1; J = 1; m_a = 1; Omega = 0:0.01:5; p_0 = 1; gamma = 0:0.25:1;

% Modulus and phase of the displacement
q_pa = zeros(length(Omega),length(gamma));
phi = zeros(length(Omega),length(gamma));

for i=1:length(gamma)
    q_pa(:,i) = 1./sqrt((1 - Omega.^2).^2 + (2*gamma(i).*Omega).^2);
    phi(:,i) = atan(2*gamma(i).*Omega./(1-Omega.^2));
end

% Plots

figure(1)
semilogy(Omega,q_pa(:,:))
axis([0 5 10E-3 50])
grid on
hold on
title("Amplitude of q_p vs dimensionless frequency","FontSize",12)
xlabel("$\Omega/\omega_0$"); ylabel("$|q_p(\Omega/\omega_0)|$")
legend

figure(2)
plot(Omega,phi(:,:))
axis([0 5 -pi/2 pi/2])
grid on
hold on
title("Phase delay of q_p vs dimensionless frequency","FontSize",12)
xlabel("$\Omega/\omega_0$"); ylabel("$\varphi (\Omega/\omega_0)$")
legend

%% TUNNED MASS DAMPER

q_tmd = zeros(length(Omega),1);
H = zeros(2,2,length(Omega));
Det = zeros(1,length(Omega));

for i  = 1:length(Omega)
    H(:,:,i) = [K+K_a-Omega(i)^2*J -K_a;
                -K_a K_a-Omega(i)^2];
    Det(:,i) = det(H(:,:,i));
end

q_tmd(:,:) = (K_a - Omega.^2*m_a)./Det;

figure(3)
plot(Omega,q_tmd(:,:))
%axis([0 5 10E-3 50])
grid on
title("Amplitude of q_p vs dimensionless frequency","FontSize",12)
xlabel("$\Omega/\omega_0$"); ylabel("$|q_p(\Omega/\omega_0)|$")
