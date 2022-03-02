%% MEF Clase 3

clc
clear
format short

%% Datos

E=1; L=1; A=1; f0=100; P = 10;

%% Datos Integración Numérica

%Cuadratura n=1
chi_ip1 = 0; %Coordenadas de los puntos de integración
w_ip1 = 2;   %Peso de integración
%Cuadratura n=2
chi_ip2 = [-1 1]*1/sqrt(3); %Coordenadas de los puntos de integración
w_ip2 = [1 1]; %Peso de integración
%Cuadratura n=3
chi_ip3 = [-sqrt(3/5) 0 sqrt(3/5)]; %Coordenadas de los puntos de integración
w_ip3 = [5/9 8/9 5/9]; %Peso de integración
%Cuadratura n=4
chi_ip4 = [-0.8611363116 -0.3399810436 0.3399810436 0.8611363116];
w_ip4 = [0.3478548451 0.6521451549 0.6521451549 0.3478548451];


ipk = 2; %Cuadratura de integración de K
if ipk==1
    chi_ipk = chi_ip1; w_ipk = w_ip1;
elseif ipk==2
    chi_ipk = chi_ip2; w_ipk = w_ip2;
elseif ipk==3
    chi_ipk = chi_ip3; w_ipk = w_ip3;
elseif ipk==4
    chi_ipk = chi_ip4; w_ipk = w_ip4;
end

ipf = 2; %Cuadratura de integración de Fuerzas
if ipk==1
    chi_ipf = chi_ip1; w_ipf = w_ip1;
elseif ipk==2
    chi_ipf = chi_ip2; w_ipf = w_ip2;
elseif ipk==3
    chi_ipf = chi_ip3; w_ipf = w_ip3;
elseif ipk==4
    chi_ipf = chi_ip4; w_ipf = w_ip4;
end


%% Matriz de rigidez

n = 2;

%Ke = Bip1'*E*A*Bip1*Je*wip1 + Bip2'*E*A*Bip2*Je*wip2
x3 = L; x1 = 0;
Le = (x3-x1);
Je = Le/2;
iJe = Je^(-1);

Ke = zeros(3);
% Be_ = zeros() %Almacén de matrices cinemáticas

for ip = 1 : ipk
    
    chi = chi_ipk(ip);
    w = w_ipk(ip);
    Be = iJe*[chi-1/2 -2*chi chi+1/2];
    Ke = Ke + Be'*E*A*Be*Je*w;
    
end
Ke;

Fe = zeros(3,1);
for ip = 1 : ipf
    
    chi = chi_ipf(ip);
    w = w_ipf(ip);
    Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
    fnod = f0*[1; 1; 1];%Valor de la fuerza aplicada en los nodos
    f_chi = Ne*fnod;
    Fe = Fe + Ne'*f_chi*Je*w;
    
end
Fe;

%% Ensamblaje
K = Ke;
Fnod = zeros(3,1);
Fnod(3) = P; %Fuerzas externas en los  nodos libres
Fneq = Fe; %
F = Fnod + Fneq;

%% Condiciones de contorno 
GDL = size(K,1);
gdl_T = [1:GDL];
gdl_R = 1;
gdl_L = setdiff(gdl_T,gdl_R);
KLL = K(gdl_L,gdl_L);
fL = F(gdl_L);
uL = KLL\fL;
U = zeros(length(gdl_T),1);
U(gdl_L) = uL;

%% Cálculo de tensiones y deformaciones
ne = 1;
nne = 3; %Número de nodos por elemento
u_e1 = U([1:3]);
epsilon_ip_e = zeros(1,ipk); %Inicialización en los elementos ...

coord_n = [0 L/2 L];
coord_n_e1 = coord_n([1 2 3]);

for i = 1:ipk
    chi = chi_ipk(i);
    Be = iJe*[chi-1/2 -2*chi chi+1/2];
    epsilon_ip_e(1,i) = Be*u_e1;
end 

sigma_ip_e = epsilon_ip_e*E; %Tensiones en los puntos de integración

% Tensiones nodales  

sigma_n_e = zeros(ne,nne); %Inicialización de las tensiones nodales en filas por elemento, en columnas por nodos
r = sqrt(3); %Razón de homotecia
chi_nip = [-1 0 1]*r;

for j=1:nne
    chi = chi_nip(j);
    Nen = [(1-chi)/2 (1+chi)/2];
    sigma_n_e(1,j) = Nen*sigma_ip_e(1,:)';
end


%% Solución fuerte 

x_ = 0:0.001:L;
for i =1:length(x_)
    x = x_(i);
    ufuerte(i) = -f0/(2*E*A)*x^2 + (P+f0*L)/(E*A)*x;
    epsfuerte(i) = -f0/(E*A)*x + (P+f0*L)/(E*A);
    sigmafuerte(i) = epsfuerte(i)*E;
    
end
    

figure(1)
plot(x_,ufuerte,'r-')

figure(2)
plot(x_,sigmafuerte,'r-')
hold on
plot(coord_n,sigma_n_e,'b*--')

coord_ip = zeros(ne,ipk); %Inicialización de las coordenadas absolutas de los puuntos de integración
for ip=1:ipk
    chi = chi_ipk(ip);
    Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
    coord_ip(1,ip) = Ne*coord_n_e1'; 
end

 plot(coord_ip,sigma_ip_e,'ks--')

% No hay continuidad en las tensiones de un elemento a otro
% Los mejores resultados en tensiones se obtienen en los puntos de
% integración
% Los mejores resultados en deformaciones se prooducen en los nodos
