%% CODE TO PRODUCE GRAPHS FOR THE TFG

clc, close, clear

%% ONE-DEGREE OF FREEDOM 

% Transfer function

K = 1; F = 1; J = 1; Omega = 1:100; p_0 = 1;

q_pa = 1./(K + 1j*Omega*F - Omega.^2*J)*p_0;
phi = atan(Omega*F./(K - Omega.^2*J));

figure(1)
plot(Omega,-q_pa)


figure(2)
plot(Omega,-phi)