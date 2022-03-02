%% Grader 2: Barra 1D con sección variable

clear
format short

%% Datos de entrada
L = 1; E = 1000; P = 10; f0 = 100; a0 = L/30;

%% Solucion analitica

x = L;
u_fuerte(x) = -3/(40*a0^2*E)*(exp(2*x/3)*(f0*(4*x^5-30*x^4+180*x^3-810*x^2+2430*x-3645)-600)+3645*f0+600);

%% Solucion MEF con elementos cuadraticos

% Datos integración numérica

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


ipk = 3; %Cuadratura de integración de K
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
if ipk==1
    chi_ipf = chi_ip1; w_ipf = w_ip1;
elseif ipk==2
    chi_ipf = chi_ip2; w_ipf = w_ip2;
elseif ipk==3
    chi_ipf = chi_ip3; w_ipf = w_ip3;
elseif ipk==4
    chi_ipf = chi_ip4; w_ipf = w_ip4;
end

At = [a0 a0*exp(-1/24) a0*exp(-1/12) a0*exp(-1/8) a0*exp(-1/6) a0*exp(-5/24) a0*exp(-1/4) a0*exp(-7/24) a0*exp(-1/3)];       % Vector de áreas total

% Rigidez local de los elementos y fuerza nodal local de los elementos

x3 = L/4; x1 = 0;
Le = x3 - x1;
Je = Le/2; iJe = 1/Je;

% Elemento 1
Ke1 = zeros(3);
for ip = 1:ipk
    chi = chi_ipk(ip);
    w = w_ipk(ip);
    Be = iJe*[chi-1/2 -2*chi chi+1/2];
    Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
    A_ip1 = Ne*At(1:3)';
    
    Ke1 = Ke1 + Be'*E*A_ip1^2*Be*Je*w;
end

% Elemento 2
Ke2 = zeros(3);
for ip = 1:ipk
    chi = chi_ipk(ip);
    w = w_ipk(ip);
    Be = iJe*[chi-1/2 -2*chi chi+1/2];
    Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
    A_ip2 = Ne*At(3:5)';   
    
    Ke2 = Ke2 + Be'*E*A_ip2^2*Be*Je*w;
end

% Elemento 3
Ke3 = zeros(3);
for ip = 1:ipk
    chi = chi_ipk(ip);
    w = w_ipk(ip);
    Be = iJe*[chi-1/2 -2*chi chi+1/2];
    Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
    A_ip3 = Ne*At(5:7)';
    
    Ke3 = Ke3 + Be'*E*A_ip3^2*Be*Je*w;
end

% Elemento 4
Ke4 = zeros(3);
for ip = 1:ipk
    chi = chi_ipk(ip);
    w = w_ipk(ip);
    Be = iJe*[chi-1/2 -2*chi chi+1/2];
    Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
    A_ip4 = Ne*At(7:9)';
    
    Ke4 = Ke4 + Be'*E*A_ip4^2*Be*Je*w;
end


% Ensamblaje: Rigidez global y fuerza nodal global

K_quest = zeros(9);

% Elemento 1 
    
gdle1 = [1 2 3];
K_quest(gdle1, gdle1) = K_quest(gdle1, gdle1) + Ke1;

% Elemento 2 

gdle2 = [3 4 5];
K_quest(gdle2, gdle2) = K_quest(gdle2, gdle2) + Ke2;

% Elemento 3 

gdle3 = [5 6 7]; 
K_quest(gdle3, gdle3) = K_quest(gdle3, gdle3) + Ke3;

% Elemento 4

gdle4 = [7 8 9];
K_quest(gdle4, gdle4) = K_quest(gdle4, gdle4) + Ke4;

% Vector de fuerzas nodales

feq_quest = zeros(9,1);

fnod = f0*[0 (1/8)^4 (1/4)^4 (3/8)^4 (1/2)^4 (5/8)^4 (3/4)^4 (7/8)^4 1];             % Vector de fuerzas en los nodos

for ip=1:ipf
   chi = chi_ipf(ip);
   w = w_ipf(ip);
   Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
   f_chi = Ne*fnod(1:3)';
   
   feq_quest(1:3) = feq_quest(1:3) + Ne'*f_chi*Je*w;
   
end

for ip=1:ipf
   chi = chi_ipf(ip);
   w = w_ipf(ip);
   Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
   f_chi = Ne*fnod(3:5)';
   
   feq_quest(3:5) = feq_quest(3:5) + Ne'*f_chi*Je*w;
   
end

for ip=1:ipf
   chi = chi_ipf(ip);
   w = w_ipf(ip);
   Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
   f_chi = Ne*fnod(5:7)';
   
   feq_quest(5:7) = feq_quest(5:7) + Ne'*f_chi*Je*w;
   
end

for ip=1:ipf
   chi = chi_ipf(ip);
   w = w_ipf(ip);
   Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
   f_chi = Ne*fnod(7:9)';
   
   feq_quest(7:9) = feq_quest(7:9) + Ne'*f_chi*Je*w;
   
end

fext = [0 0 0 0 0 0 0 0 P];         % Fuerzas externas

ft = feq_quest + fext';              % Vector de fuerzas totales

% Condiciones de contorno: calcula desplazamientos nodales (Condensación estática)

GDL = size(K_quest,1);
gdl_T = [1:GDL];
gdl_L = [2:9]; gdl_R = [1];
KLL = K_quest(gdl_L, gdl_L);
KRR = K_quest(gdl_R,gdl_R);
KRL = K_quest(gdl_R,gdl_L);
KLR = KRL';

fL = ft(gdl_L,1);
uL = KLL\fL;
uL_quest = uL(8);                   % Desplazamiento en nodos libres (el extremo libre)
U = zeros(length(gdl_T),1);
U(gdl_L) = uL;

% El error relativo % en desplazamientos en el extremo libre entre la solución fuerte y la solución mef (error_u)

error_abs = abs(u_fuerte-uL(8));
error_u = error_abs/u_fuerte*100;

%% Determinar la deformacion y tension

ne = 4;                    % número de elementos 
nne = 3;                   % número de nodos por elemento
r = sqrt(3);               % razón de homotecia
chi_nip = [-1 0 1]*r;

coord_n = [0 L/8 L/4 3*L/8 L/2 5*L/8 3*L/4 7*L/8 L];                    % coordenadas de los nodos totales
coord_n_e1 = coord_n([1 2 3]);                                          % coordenadas de los nodos de cada elemento
coord_n_e2 = coord_n([3 4 5]);  
coord_n_e3 = coord_n([5 6 7]);
coord_n_e4 = coord_n([7 8 9]);

ipk_ = 2;                                                               % Redefino el orden porque tengo 2 puntos de integración chi

% Elemento 1

u_e1 = U([1:3]);
epsilon_ip_e1 = zeros(1,ipk_);

for i = 1:ipk_
    chi = chi_ip2(i);
    Ben = iJe*[chi-1/2 -2*chi chi+1/2];
    epsilon_ip_e1(1,i) = Ben*u_e1;
end 

sigma_ip_e1 = epsilon_ip_e1*E; % Tensiones en los puntos de integración

sigma_n = zeros(ne,nne);       % Inicialización de las tensiones nodales en filas por elemento, en columnas por nodos

for j=1:nne
    chi = chi_nip(j);
    Nen =[(1-chi)/2 (1+chi)/2];                       
    sigma_n(1,j) = Nen*sigma_ip_e1(1,:)';
end

% Elemento 2

u_e2 = U([3:5]);
epsilon_ip_e2 = zeros(1,ipk_);

for i = 1:ipk_
    chi = chi_ip2(i);
    Ben = iJe*[chi-1/2 -2*chi chi+1/2];
    epsilon_ip_e2(1,i) = Ben*u_e2;
end 

sigma_ip_e2 = epsilon_ip_e2*E;

for j=1:nne
    chi = chi_nip(j);
    Nen = [(1-chi)/2 (1+chi)/2];                        
    sigma_n(2,j) = Nen*sigma_ip_e2(1,:)';
end

% Elemento 3

u_e3 = U([5:7]);
epsilon_ip_e3 = zeros(1,ipk_);

for i = 1:ipk_
    chi = chi_ip2(i);
    Ben = iJe*[chi-1/2 -2*chi chi+1/2];
    epsilon_ip_e3(1,i) = Ben*u_e3;
end 

sigma_ip_e3 = epsilon_ip_e3*E; 

for j=1:nne
    chi = chi_nip(j);
    Nen = [(1-chi)/2 (1+chi)/2];                        
    sigma_n(3,j) = Nen*sigma_ip_e3(1,:)';
end

% Elemento 4

u_e4 = U([7:9]);
epsilon_ip_e4 = zeros(1,ipk_);

for i = 1:ipk_
    chi = chi_ip2(i);
    Ben = iJe*[chi-1/2 -2*chi chi+1/2];
    epsilon_ip_e4(1,i) = Ben*u_e4;
end 

sigma_ip_e4 = epsilon_ip_e4*E; 

for j=1:nne
    chi = chi_nip(j);
    Nen = [(1-chi)/2 (1+chi)/2];                        
    sigma_n(4,j) = Nen*sigma_ip_e4(1,:)';
end

epsilon_ip_t = [epsilon_ip_e1; epsilon_ip_e2; epsilon_ip_e3; epsilon_ip_e4];        % Matriz con las deformaciones en los puntos de integración
sigma_ip_t = [sigma_ip_e1; sigma_ip_e2; sigma_ip_e3; sigma_ip_e4];                  % Matriz con las tensiones en los puntos de integración

%% Solución fuerte para dibujar

x_ = 0:0.01:L;
for i =1:length(x_)
    x = x_(i);
    ufuerte(i) = -3/(40*a0^2*E)*(exp(2*x/3)*(f0*(4*x^5-30*x^4+180*x^3-810*x^2+2430*x-3645)-600)+3645*f0+600);
    epsfuerte(i) = -((f0*x^5-150)*exp((2*x)/3))/(5*E*a0^2);
    sigmafuerte(i) = epsfuerte(i)*E;
    
end

%% Errores

% Solución fuerte para el error

x_e = 0:0.125:L;
for i =1:length(x_e)
    x = x_e(i);
    ufuerte_error(i) = -3/(40*a0^2*E)*(exp(2*x/3)*(f0*(4*x^5-30*x^4+180*x^3-810*x^2+2430*x-3645)-600)+3645*f0+600);
    epsfuerte_error(i) = -((f0*x^5-150)*exp((2*x)/3))/(5*E*a0^2);
    sigmafuerte_error(i) = epsfuerte_error(i)*E;
    
end

sf_max = max(sigmafuerte_error);
sn_max = max(sigma_n(:));
sip_max = max(sigma_ip_t(:));

% Error de la solución MEF en los nodos

err_sigmamax_nod = abs(sf_max - sn_max)/sf_max*100;


% Error de la solución MEF en los puntos de integración

err_sigmamax_ip = - abs(max(sigmafuerte) - sip_max)/max(sigmafuerte)*100;

%% Representación gráfica

figure(1)
plot(x_,ufuerte,'r-')

figure(2)
plot(x_,sigmafuerte,'r-')
hold on
plot(coord_n(1:3),sigma_n(1,:),'b*--')
plot(coord_n(3:5),sigma_n(2,:),'b*--')
plot(coord_n(5:7),sigma_n(3,:),'b*--')
plot(coord_n(7:9),sigma_n(4,:),'b*--')

figure(3)
grafstress_strong = plot(x_,sigmafuerte,'r-','LineWidth',2);

figure(4)
grafstressnod = plot(coord_n(1:3),sigma_n(1,:),'go-','LineWidth',2); %representa el elemento 1 y almacena en la variable grafstressnod
hold on %mantiene los datos de la gráfica anterior cuando representes los siguientes elementos
grafstressnod = plot(coord_n(3:5),sigma_n(2,:),'go-','LineWidth',2); %representa el elemento 2
grafstressnod = plot(coord_n(5:7),sigma_n(3,:),'go-','LineWidth',2); %representa el elemento 3
grafstressnod = plot(coord_n(7:9),sigma_n(4,:),'go-','LineWidth',2); %representa el elemento 4

title('Tension en solución fuerte vs MEF en nodos','Fontsize',10);


% coord_ip = zeros(ne,ipk); %Inicialización de las coordenadas absolutas de los puuntos de integración
% for ip=1:ipk_
%     chi = chi_ip2(ip);
%     Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
%     coord_ip(1,ip) = Ne*coord_n_e1'; 
% end
% 
% plot(coord_ip,sigma_ip_t,'ks--')


% function A = A(x, a0,L)
%     A = (a0*exp(-x/(3*L)))^2
% end