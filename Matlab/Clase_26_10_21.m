%% Clase 26/10/21
clear
clc
%% Problema de la viga hermítica de 2 nodos (grados flectores)

E = 1; A = 100; I = 1; L = 1; f0y = -1; Py = 10;

%% Parámetros de integración
% cuadratura n=1
chi_ip1 = 0; % coordenadas de los puntos de integración
w_ip1=2; % peso de integración

% cuadratura n=2
chi_ip2 = [-1 1]*1/sqrt(3); % coordenadas de los puntos de integración
w_ip2 = [1 1]; % peso de integración

% cuadratura n=3
chi_ip3 = [-sqrt(3/5) 0 sqrt(3/5)]; % coordenadas de los puntos de integración
w_ip3 = [5/9 8/9 5/9]; % peso de integración

% cuadratura n=4
chi_ip4 = [-0.8611363116 -0.3399810436 0.3399810436 0.8611363116];
w_ip4 = [0.3478548451 0.6521451549 0.6521451549 0.3478548451];

ipK = 2; % cuadratura de integración de K
if ipK == 1
    chi_ipK = chi_ip1;
    w_ipK = w_ip1;
elseif ipK == 2
    chi_ipK = chi_ip2;
    w_ipK = w_ip2;
elseif ipK == 3
    chi_ipK = chi_ip3;
    w_ipK = w_ip3;
elseif ipK == 4
    chi_ipK = chi_ip4;
    w_ipK = w_ip4;
end

ipf = 3; % cuadratura de integración de K
if ipf == 1
    chi_ipf = chi_ip1;
    w_ipf = w_ip1;
elseif ipf == 2
    chi_ipf = chi_ip2;
    w_ipf = w_ip2;
elseif ipf == 3
    chi_ipf = chi_ip3;
    w_ipf = w_ip3;
elseif ipf == 4
    chi_ipf = chi_ip4;
    w_ipf = w_ip4;
end


%% Matriz de rigidez
n = 2;

x2 = L;
x1 = 0;
Le = x2-x1;
Je = Le/2;
iJe = 1/Je; % inversa del jacobiano

Kef = zeros(4); % inicializacion de la matriz de rigidez

 for ip = 1:ipK
     chi = chi_ipK(ip);
     w = w_ipK(ip);
     k1 = (chi*3)/2; k2 = (chi*3-1)/2; k3 = -(chi*3)/2; k4 = (chi*3+1)/2;
     Bef = iJe^2*[k1 Je*k2 k3 Je*k4]; % matriz cinemática de flexión
     
     Kef = Kef +  Bef'*E*I*Bef*Je*w; % Matriz de rigidez de flexión
 end

Kef

% si uso dos nodos para las fuerzas tengo que usar las funciiones de forma
% de grado mínimo necesario (lineales si usamos dos nodos)


%% Vectores de fuerzas nodales equivalentes
Feqf = zeros(4,1); % inicialización del vector de fuerzas nodales 
% 4 pq hay 4 gdl en la viga de 2 nodos

for ip = 1:ipf
     chi = chi_ipf(ip);
     w = w_ipf(ip);
     
     nh1 = -0.75*chi + 0.25*chi^3 + 0.5;
     nh2 = -0.25*chi - 0.25*chi^2 + 0.25*chi^3 + 0.25;
     nh3 = 0.75*chi - 0.25*chi^3 + 0.5;
     nh4 = -0.25*chi + 0.25*chi^2 + 0.25*chi^3 - 0.25;
     
     Nhe = [nh1 Je*nh2 nh3 Je*nh4]; % Funciones hermíticas
     
     fy0_ = [0;f0y];
     NLe = [(1-chi)/2 (1+chi)/2];
     
     f_chi = NLe*fy0_;
     
     Feqf = Feqf + Nhe'*f_chi*Je*w;
end
 
Feqf













