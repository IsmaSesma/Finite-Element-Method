%% Grader 3: Viga hermítica de sección variable con 4 elementos

clc;clear;close all;
format short
%% Datos de entrada

E = 210E9;  nu = 0.3; L = 0.7; a0 = 0.016; t = a0; q0 = 0; P = -50;

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

ipf = 3; %Cuadratura de integración de Fuerzas
if ipf==1
    chi_ipf = chi_ip1; w_ipf = w_ip1;
elseif ipf==2
    chi_ipf = chi_ip2; w_ipf = w_ip2;
elseif ipf==3
    chi_ipf = chi_ip3; w_ipf = w_ip3;
elseif ipf==4
    chi_ipf = chi_ip4; w_ipf = w_ip4;
end

%% Matriz de rigidez de la viga

ne = 10;   % número de elementos
nn = 11;   % número de nodos
gdln = 3;  %grados de libertad por nodo
GDL = 33;  %grados de libertad totales

coord_n = zeros(nn,2); % Matriz de coordenadas nodales: en filas nodos, en columna 1 coordenada x, en columna 2 coordenada y
                       % Numero de filas = numero de nodos, 2 columnas
coord_n(:,1) = 0:L/ne:L;    % Solo hay coordenada x en los nodos

conectividad_e = zeros(ne,2); % Matriz de conectividad de elementos a través de los nodos
conectividad_e(:,1) = 1:1:ne;
conectividad_e(:,2) = 2:1:nn;

% Bucle de integracion y montaje de elementos

K = zeros(nn*gdln);   % Declara la matriz de rigidez global
F = zeros(nn*gdln,1); % Declara el vector de fuerzas nodales en c.globales

for e = 1:ne % Almacen de los grados de libertad de cada elemento   
    index = conectividad_e(e,:); % Indice de posicion de los nodos en global
    x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
    Le = x2-x1; Je = Le/2; iJe = 1/Je; % Jacobiano de la transformacion de cada elemento
    gdle = [index(1)*gdln-2 index(1)*gdln-1 index(1)*gdln...
            index(2)*gdln-2 index(2)*gdln-1 index(2)*gdln]; %indince de cada elemento
    gdlea =  gdle([1,4]);     % Grados de libertad axiles
    gdlef =  gdle([2,3,5,6]); % Grados de libertad flectores
    
    % Matriz de rigidez axil
    
    Kea = zeros(2); 
    for i = 1:1 % Integracion numerica de Kea
        % Se hace como en los graders de barras
        chi = chi_ip1(i); w = w_ip1(i);
        Be = iJe*[-1/2 1/2]; % Para 2 nodos
        Ne = [(1-chi)/2 (1+chi)/2]; % Para 2 nodos
        
        Ae_ = [A_(x1,a0,L,t); A_(x2,a0,L,t)];
        Ae = Ne*Ae_;
        
        Kea = Kea + Be'*E*Ae*Be*Je*w;
    end
    K(gdlea,gdlea)= K(gdlea,gdlea)+Kea; % Ensamblaje de la rigidez axil en la matriz global
    
    % Matriz de rigidez a flexion
    
    Kef = zeros(4);
    
    for i = 1:2 % Integracion numerica de Kef
        % Como lo hicimos en clase
        chi = chi_ip2(i); w = w_ip2(i);
        Ne = [(1-chi)/2 (1+chi)/2];             % Tengo 2 nodos
        k1 = (chi*3/2); k2 = (chi*3-1)/2 ; k3 = -(chi*3/2); k4 = (chi*3+1)/2;
        Bef = iJe^2*[k1 Je*k2 k3 Je*k4];        % Matriz cinematica de flexion
        
        Ie_ = [Iy_(x1,a0,L,t);Iy_(x2,a0,L,t)];
        Ie = Ne*Ie_;
        
        Kef = Kef + Bef'*E*Ie*Bef*Je*w; 
        
    end
    
    K(gdlef,gdlef)=K(gdlef,gdlef)+Kef; % Ensamblaje de la rigidez de flexion en la matriz global
    
    
    % Fuerzas axiles

    fea = zeros(2,1);
    
    for i = 1:2 % Integracion numerica de fea
        % Como en  Grader 2
        chi = chi_ip2(i); w = w_ip2(i);
        Ne = [(1-chi)/2 (1+chi)/2]; % Para 2 nodos
        fea_ = [0;0];
        fea = fea+Ne'*Ne*fea_*Je*w;
    end
     
    F(gdlea)= F(gdlea)+fea; % Ensamblaje de la rigidez axil en la matriz global
     
    % Fuerzas de flexion
     
    fef = zeros(4,1);
     
    for i = 1:3 % Integracion numerica de fef
        % Como en clase de flexión
        chi = chi_ip3(i); w = w_ip3(i);
        NLe = [(1-chi)/2 (1+chi)/2]; % Para 2 nodos;
        nh1 = -0.75*chi+0.25*chi^3+0.5; nh2 = -0.25*chi-0.25*chi^2+0.25*chi^3+0.25;
        nh3 = 0.75*chi-0.25*chi^3+0.5; nh4 = -0.25*chi+0.25*chi^2+0.25*chi^3-0.25;
        Nhe = [nh1 Je*nh2 nh3 Je*nh4]; % Funciones hermiticas del elemento (de forma)
        
        fef_= [0;0]; 
        
        fef = fef+Nhe'*NLe*fef_*Je*w;
    end
    F(gdlef) = F(gdlef)+fef; % Ensamblaje de la rigidez de flexion en la matriz global
end

K_quest = K;

%% Cálculo de u1_quest y sigma1_quest

p = zeros(nn*gdln,1);        % Cargas externas aplicadas en los nodos
p(33) = P;
ft = F + p;                  % Fuerzas totales   

u = zeros(nn*gdln,1);
gdl_L = [4:nn*gdln]; gdl_R = [1 2 3]; % Grados libres (gdl_L) y restringidos (gdl_R)

K_LL = K(gdl_L,gdl_L);
K_LR = K(gdl_L,gdl_R);

f_L = ft(gdl_L); % Cargas en los nodos libres
u_L = K_LL\f_L;

u(gdl_L) = u_L; % Vector de desplazamientos total

f_R = K_LR'*u_L; % Reacciones

u1_quest = u(13:15);

%% Cálculo u2_quest

Baric = 3*L/8;                   % Baricentro de la distribución       
Coord_P = [Baric 0];             % Lo uso para el bucle

for i=1:ne 
    index = conectividad_e(i,:); % Índice de posición de los nodos en global 
    x1 = coord_n(index(1),1); x2 = coord_n(index(2),1); 
    gdle = [index(1)*gdln-2 index(1)*gdln-1 index(1)*gdln...                                
        index(2)*gdln-2 index(2)*gdln-1 index(2)*gdln]; 
    gdlea =  gdle([1,4]); gdlef =  gdle([2,3,5,6]); 
    if Coord_P(1)>x2
        continue 
    elseif or(Coord_P(1)==x2,Coord_P(1)==x1)     
        disp("Los desplazamientos ya están calculados")
        break 
    elseif Coord_P(1)<x2
        break 
    end
end

%Coordenada xi y datos para el cálculo

xi_p = 2/(x2 - x1)*(Baric - (x1+x2)/2);
Le = x2 - x1; Je = Le/2; iJe = 1/Je;
NH = Nh(xi_p,Je);
NL = N_2_nod(xi_p);
dNH = dNh(xi_p,Je);

Bf = Bh(xi_p,Je); % Matriz cinemática de flexión

% Desplazamientos (axil, flecha y giro) =====> Diapo 42

ua_p = NL*u(gdlea);
uf_p = NH*u(gdlef);
ug_p = iJe*dNH*u(gdlef);

% Vector Desplazamiento

u_p = [ua_p ;uf_p ;ug_p];
u2_quest = u_p;

%% Cálculo de tensiones 

for e = 1:ne
    
    index = conectividad_e(e,:); % Indice de posicion de los nodos en global
    x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
    u_e(:,e) = [u([2+3*(e-1); 3+3*(e-1); 5+3*(e-1); 6+3*(e-1)])]; % Cojo los DOF salvo los de axil    

% Deformaciones y tensiones en los puntos de integración (Cómo en Grader2 pero para flexión)

    epsilon_ip_e = zeros(ipk,1); %deformación del elemento en los puntos de integración de K
    
    for ip = 1:2
        chi = chi_ip2(ip);  % Datos de integración    
        k1=(chi*3/2); k2=(chi*3-1)/2; k3=-(chi*3/2); k4=(chi*3+1)/2;
        Bef = iJe^2*[k1 Je*k2 k3 Je*k4]; % Matriz cinemática
        Ne = [(1-chi)/2 (1+chi)/2];
        ae_= [a_(x1,a0,L); a_(x2,a0,L)];
        z_i = 0.5*Ne*ae_;   % z es 1/2*a(x)
        epsilon_ip_e(ip) = -z_i*Bef*u_e(:,e);    
    end
    
    sigma_ip(:,e) = epsilon_ip_e*E;
    
% Tensiones en los nodos
    
    chi_1 = -sqrt(3); % Coordenadas locales de los nodos
    chi_2 = sqrt(3);  % sqrt(3) porque se usan dos chis por elemento
    
    sigma_n(:,e) = [[(1-chi_1)/2 (1+chi_1)/2]*sigma_ip(:,e), [(1-chi_2)/2 (1+chi_2)/2]*sigma_ip(:,e)];
   
end

sigma1_quest = max(max(sigma_n))/1E6; % En MPa

chi_baric = 0;
sigma2_quest = - [(1-chi_baric)/2 (1+chi_baric)/2]*sigma_n(:,2)/1E6;


%% Errores

% Necesito la solución analítica ==> Ecuación de Navier

x_ = 0:0.001:L;
for i=1:length(x_)
    x = x_(i);
    My(i) = (0.5*q0*x^2*(1-(1/6)*(x/L)^2)+(-P-(2*q0*L/3))*x + (q0*L^2/4)+P*L);
    a(i) = a0*(1-(2*x)/(3*L));
    b(i) = 3*a(i);
    Iy (i)= ((b(i)-t)*t^3)/6 + 2*(b(i)-t)*t*((a(i)-t)/2)^2 +(t*(a(i)-t)^3)/6;
    z(i) = a(i)/2;
    sigmafuerte(i) =  (My(i)*z(i)/Iy(i));
end

s1_max = max(abs(sigmafuerte(:)))/1E6;
s2_max = max(abs(sigmafuerte(375)))/1E6;

err1_quest = abs(s1_max - sigma1_quest)/s1_max*100;
err2_quest = abs(s2_max + sigma2_quest)/s2_max*100;

%%

sigma_ip = zeros(2,ne);
x_ip = zeros(2,ne);
%
for e = 1:ne; % Indice de posicion de los nodos en global  
    index = conectividad_e(e,:); % Indice de posicion de los nodos en global
    x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
    gdle = [index(1)*gdln-2 index(1)*gdln-1 index(1)*gdln...
            index(2)*gdln-2 index(2)*gdln-1 index(2)*gdln];
    gdlea =  gdle([1,4]); % Grados axiles
    gdlef =  gdle([2,3,5,6]); % Grados flectores
    %
    Le = x2-x1; Je = Le/2; iJe = 1/Je; % Jacobiano de la transformacion de cada elemento
    %
    for ip = 1:2
        %
        chi = chi_ip2(ip);
        %
        NL = N_2_nod(chi);
        NH = Nh(chi,Je);
        d_NH = dNh(chi,Je);
        %
        x_ip(ip,e) = NL*[x1;x2];
        %
        ua_ip(ip,e) = NL*u(gdlea);
        uf_ip(ip,e) = NH*u(gdlef);
        ug_ip(ip,e) = iJe*d_NH*u(gdlef);
        %
        % Matriz 3D para guardar desplazamientos y giros de los ptos de integracion
        % (1,2,3)   1 -> Guarda valores de un pto de integracion
        %           2 -> Pto de integracion
        %           3 -> Elemento
        %
        u_ip(:,ip,e) = [ua_ip(ip,e); uf_ip(ip,e); ug_ip(ip,e)];
        Bf = Bh(chi,Je);
        K_ip(ip,e) = Bf*u(gdlef);
        I_ip(ip,e) = NL*[Iy_(x1,a0,L,t);Iy_(x2,a0,L,t)];
        M_ip(ip,e) = E*I_ip(ip,e)*K_ip(ip,e);
        zmax_ip(ip,e) = NL*[zmax_(x1,a0,L);zmax_(x2,a0,L)];
        sigma_ip(ip,e) = E*K_ip(ip,e)*zmax_ip(ip,e);
        %
    end
end
%

sigma_nod = zeros(2,ne);
chi_nodal = [-1;1];
%
for e = 1:ne
    index = conectividad_e(e,:); % Indice de posicion de los nodos en global
    x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
    gdle = [index(1)*gdln-2 index(1)*gdln-1 index(1)*gdln...
            index(2)*gdln-2 index(2)*gdln-1 index(2)*gdln];
    gdlea =  gdle([1,4]); % Grados axiles
    gdlef =  gdle([2,3,5,6]); % Grados flectores
    %
    chi_1 = chi_ip2(1); chi_2 = chi_ip2(2);
    %
    N1 = N_2_nod(chi_1); N2 = N_2_nod(chi_2);
    %
    M_Interpolacion = [N1;N2];
    %
    sigma_nod(:,e) = M_Interpolacion\sigma_ip(:,e);
    M_nod(:,e) = M_Interpolacion\M_ip(:,e);
    %
end

x_L=0:L/10:L;
sigmaf_n = zeros(ne,2) %matriz de almacen de sigma en los nodos, en filas por elementos y en columnas por nodos
Mze_n = zeros(ne,2) %matriz de almacen de sigma en los nodos, en filas por elementos y en columnas por nodos
Mze_ = M_nod'

for e = 1:ne
    
x_e= zeros(1,2);
x_e=x_L(e:e+1);
    


%calcula la tension en los nodos de cada elemento y elemento a elemento
sigmaf_n(e,:) = sigma_ip(:,e)'/(10^6)
Mze_n(e,:) = M_ip(:,e)'

figure(1)
graf1=plot(x_e,Mze_n(e,:),'ro--','LineWidth',1); % x_n son las coordenadas de los nodos, y tau_n es la tensión en MPa
hold on
title('Momentos MEF en nodos','Fontsize',10);

figure(2)
graf2=plot(x_e,sigmaf_n(e,:),'ro--','LineWidth',1);   % x_n son las coordenadas de los nodos, y sigmaf_n es la tensión en MPa
hold on
title('Tension MEF en nodos','Fontsize',10);

end

%% Representa desplazamientos
Ux = u(1:3:(3*nn));
Uy = u(2:3:(3*nn)); Uymax = max(abs(Uy)); %flecha maxima
Uxy = u(3:3:(3*nn));
%escala = 0.15*L/Uymax;
coordx_n=coord_n(:,1);
coordy_n=coord_n(:,2);
figure(3)
graf3 = plot(coordx_n,coordy_n,'bo-',coordx_n+Ux,coordy_n+Uy,'rs--','LineWidth',2)
axis([0 (L+max(Ux)*10) min(Uy)*2 max(Uy)*2])
hold on
title("Desplazamientos nodales",'Fontsize',12)
xlabel("coordenada nodal X"); ylabel("coordenada nodal Y");


%% Funciones

function Iy = Iy_(x,a0,L,t)  % Inercia de la seccion
    a = a0*((-2/(3*L))*x+1); % Lado variable
    b = 3*a;                 % Ancho variable
    Iy = (1/12)*b*a^3-(1/12)*(b-2*t)*(a-2*t)^3;
end

% Altura en los nodos

function a = a_(x,a0,L)  
    a = a0*((-2/(3*L))*x+1);
end

% Cálculo del Area

function A = A_(x,a0,L,t) % Area de la seccion
    a = a0*exp(-x/L); % Canto variable
    b = 3*a; % Ancho variable
    A = a*b-(a-2*t)*(b-2*t); % Area de la seccion
%     A = 2*t*(a+b-2*t); % Area de la seccion    
end

% Funciones de forma

function Ne = N_2_nod(chi) %Funcion de forma de orden 2
Ne = [(1-chi)/2 (1+chi)/2];
end

% Funciones de forma hermíticas

function NH = Nh(xi,Je)    
h1= (2-3*xi+xi^3)*1/4;    
h2= (1-xi-xi^2+xi^3)*1/4;    
h3= (2+3*xi-xi^3)*1/4;    
h4= (-1-xi+xi^2+xi^3)*1/4;
    NH = [h1 Je*h2 h3 Je*h4];
end

% Derivada primera de la función de forma

function NH_prima = dNh(xi,Je)
h1= (-3+3*xi^2)*1/4;    
h2= (-1-2*xi+3*xi^2)*1/4;    
h3= (3-3*xi^2)*1/4;    
h4= (-1+2*xi+3*xi^2)*1/4;    
NH_prima = [h1 Je*h2 h3 Je*h4];
end

% Matriz cinemática de las funciones de forma

function BH = Bh(xi,Je)    
iJe = Je^-1;    
k1 =(xi*3/2); k2=(xi*3/2-1/2); k3 = (-xi*3/2); k4=(xi*3/2+1/2);    
BH = iJe^2*[k1 Je*k2 k3 Je*k4];
end

function Be = B_3_nod(chi,Je)    
Be = (1/Je)*[chi-1/2 -2*chi chi+1/2];    
end 
       
% Carga distribuida

function q = qp(x,q0,L)
    q = q0*(1-(x/L)^2);
end

function z = zmax_(x,a0,L)
a = a0*((-2/(3*L))*x+1);   
z = a/2;
end      
        