%% Integración numérica para elemento barra de 3 nodos

%% Datos

E = 1;
A = 1;
L = 1;
f0 = 1;

%% Parámetros de integración

% n=1 
chi_ip1 = 0; % coordenadas de los puntos de integración
w_ip1=2;     % peso de integración

% n=2
chi_ip2 = [-1 1]*1/sqrt(3); % coordenadas de los puntos de integración
w_ip2 = [1 1];              % peso de integración

% n=3
chi_ip3 = [-sqrt(3/5) 0 sqrt(3/5)]; % coordenadas de los puntos de integración
w_ip3 = [5/9 8/9 5/9];              % peso de integración


%% Matriz de rigidez

x3 = L;
x1 = 0;
Le = x3-x1;
Je = Le/2;
iJe = 1/Je;     % inversa del jacobiano

Ke = zeros(3);  % inicializacion de la matriz de rigidez

 for ip = 1:2
     chi = chi_ip2(ip);
     w = w_ip2(ip);
     Be = iJe*[chi-1/2 -2*chi chi+1/2];
     
     Ke = Ke +  Be'*E*A*Be*Je*w;
 end

Ke


%% Vectores de fuerzas nodales equivalentes

Fe = zeros(3,1); % inicialización del vector de fuerzas nodales 

for ip = 1:2
     chi = chi_ip2(ip);
     w = w_ip2(ip);
     Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
     fnod = f0*[1;1;1];
     f_chi = Ne*fnod;
     
     Fe = Fe + Ne'*f_chi*Je*w;
end
 
Fe




