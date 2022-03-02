%% Grader 4: Elementos viga Timoshenko y MPCs

clear
format long g
%% Datos

% Geométricos generales

E = 70E9; L = 5; nu = 0.3; G = E/(2*(1 + nu));
u = - L/100;    % Desplazamiento dato en el nodo 4 (hacia abajo)

% Para la integracion numerica (lo de siempre) 

chi_ip1 = 0; w_ip1 = 2; % n=1
chi_ip2 = [-1 1]*1/sqrt(3); w_ip2 = [1 1]; % n=2
chi_ip3 = [-sqrt(3/5) 0 sqrt(3/5)]; w_ip3 = [5/9 8/9 5/9]; % n=3


% Para cálculo en los bucles

nne = 2;            % Número de nodos de cada elemento
nn_v1 = 17;         % Nodos de la primera viga
nn_r = 5;           % Nodos de la riostra
ne_v1 = 16;         % Número de elementos en la viga 1
ne_r = 4;           % Número de elementos de la riostra

gdln_v1 = 3;        % Grados de liberad por nodo en la viga 1
gdln_r = 3;         % Grados de libertad por nodo en la riostra
GDL_v1 = nn_v1*gdln_v1;  % Grados de libertad de la viga 1
GDL_r = nn_r*gdln_r;     % Grados de liberad de la riostra
GDL = GDL_v1 + GDL_r;    % Grados de libertad de todo el problema

% Datos para los bucles (Ver grader 3)

coord_vl = zeros(nn_v1,2); % Coordenadas de los nodos de la viga 1
coord_vl(:,1) = 0:L/ne_v1:L; % Solo hay coordenada x en los nodos

coord_r = zeros(nn_r,2); % Coordenadas de los nodos de la riostra
% Aquí la viga está girada y tengo diferentes coordenadas, tanto x como y
coord_r(:,1) = 0:(L/3)/ne_r:L/3;
coord_r(:,2) = -L/2:(L/2)/ne_r:0;

conectividad_e = zeros(ne_v1,2); % Matriz de conectividad de elementos de la viga a través de los nodos
conectividad_e(:,1) = 1:1:ne_v1;
conectividad_e(:,2) = 2:1:nn_v1;

% Elemento viga 1

a1 = L/40; b1 = a1/2; t1 = a1/10;
kc1 = 2.5;      % Coeficiente de cortante en la viga 1
I1 = 1/12*t1*(a1 - 2*t1)^3 + 2*((1/12)*b1*t1^3 + b1*t1*((a1 - t1)/2)^2);  % Momento de inercia de la viga 1
A1 = 2*t1*b1 + t1*(a1 - 2*t1);          % Área de la sección de la viga 1

% Elemento viga 2 (la riostra)

r2 = a1/5; t2 = r2/5;
kc2 = 2;        % Coeficiente de cortante en la riostra
I2 = pi/4*(r2^4 - (r2 - t2)^4);         % Momento de inercia de la riostra
A2 = pi*(r2^2 - (r2 - t2)^2);           % Área de la sección de la riostra

%% Integración y ensamblaje (misma idea para bucles que en grader 3)

K = zeros(GDL);             % Matriz de rigidez global
F = zeros(GDL,1);           % Vector de fuerzas nodales en coordenadas globales

% Para viga 1 (copiar y pegar del grader anterior)

for e = 1:ne_v1
    coord_n = coord_vl;     % Lo renombro porque estoy usando la variable del grader 3
    gdln = gdln_v1;
    
    index = conectividad_e(e,:); % Indice de posicion de los nodos en global
    x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
    Le = x2-x1; Je = Le/2; iJe = 1/Je; % Jacobiano de la transformacion de cada elemento
    gdle = [index(1)*gdln-2 index(1)*gdln-1 index(1)*gdln...
            index(2)*gdln-2 index(2)*gdln-1 index(2)*gdln]; %indince de cada elemento
    gdlea =  gdle([1,4]);     % Grados de libertad axiles
    gdlef =  gdle([2,3,5,6]); % Grados de libertad flectores
    
    % Matriz de rigidez axil
    
    Kea_v1 = zeros(2);
    
    for i = 1:1
        % Como en los grader de barras
        chi = chi_ip1(i); w = w_ip1(i);
        Be = iJe*[-1/2 1/2]; % Para 2 nodos
        Ne = [(1-chi)/2 (1+chi)/2]; % Para 2 nodos
        
        Ae_ = [A1; A1];
        Ae = Ne*Ae_;
        
        Kea_v1 = Kea_v1 + Be'*E*Ae*Be*Je*w;
    end
    
    K(gdlea,gdlea)= K(gdlea,gdlea) + Kea_v1;
    
end 

Aef_v1 = A1/kc1;        % Área efectiva a cortadura (lo hacemos con coeficiente en vez de calcularlo
                        % como en estructuras)

Kef_v1 = Kbeam(Le,I1,Aef_v1,E,G);

% Ensamlaje

Kef = zeros(GDL);

for i = 2:3:47;
    Kef([i i+1 i+3 (i+3)+1],[i i+1 i+3 (i+3)+1]) = Kef([i i+1 i+3 (i+3)+1],[i i+1 i+3 (i+3)+1]) + Kef_v1;
end

% Como Timoshenko no considera el axil lo añado
K = K + Kef;           % Añado a la rigidez axil la flexión/cortadura por Timoshenko

% Para la riostra (la considero también una viga)

Le_r = sqrt((L/3)^2+(L/2)^2)/ne_r;          % Longitud de cada elemento de la riostra
Je_r = Le_r/2;                              % Jacobiano de la transformación
iJe = Je_r^(-1);                            % Inversa del jacobiano

% Matriz de rigidez axil de la riostra

Aef_r = A2/kc2;                             % Área efectiva de la riostra
Kef_r = Kbeam(Le_r, I2, Aef_r, E, G);
Kea_r = zeros(nne*gdln_r);

for i = 1:1
    
    Je = Je_r;
    chi = chi_ip1(i); w = w_ip1(i);
    Be = iJe*[-1/2 1/2]; % Para 2 nodos
    Ne = [(1-chi)/2 (1+chi)/2]; % Para 2 nodos
	
    Ae_ = [A2; A2];
    Ae = Ne*Ae_;
	
    Kea_r([1 4],[1 4]) = Kea_r([1 4],[1 4]) + Be'*E*Ae*Be*Je*w; % [1,4] porque son los gdl axiles

end
    
% Repito el procemiento de añadir a Timoshenko el axil
Ke_r = zeros(nne*gdln_r);       % Matriz de rigidez de un elemento de la riostra
Ke_r([2 3 5 6],[2 3 5 6]) = Ke_r([2 3 5 6],[2 3 5 6]) + Kef_r;
Ke_riostra = Kea_r + Ke_r;          % No me vale porque tengo que girar

% Cambio de base

alfa = atan((L/2)/(L/3));
R_ = [cos(alfa) sin(alfa) 0; -sin(alfa) cos(alfa) 0; 0 0 1];
R = [R_ zeros(3); zeros(3) R_];
Ke_r_giro = R'*Ke_riostra*R;    % Matriz de rigidez de la riostra girada

% Ensamblaje

Ke_r_f = zeros(nn_r*gdln_r);

for i = 1:3:10
    Ke_r_f([i i+1 i+2 i+3 i+4 i+5],[i i+1 i+2 i+3 i+4 i+5]) = Ke_r_f([i i+1 i+2 i+3 i+4 i+5],[i i+1 i+2 i+3 i+4 i+5]) + Ke_r_giro;
end
% Ensamblaje de la global

K([GDL_v1+1:GDL],[GDL_v1+1:GDL]) = K([GDL_v1+1:GDL],[GDL_v1+1:GDL]) + Ke_r_f; % Sumo a la K las componentes de la riostra, de ahí el + 1

K_quest = K; 

%% Matriz de coeficientes de restricción

p = 10;             % Restricciones = 9 por apoyos + 1 del desplazamiento u
C = zeros(p,size(K,1));

% SPCs (como explicado en clase y en ejemplo)

% Apoyo simple en 1
C(1,1) = 1; % u1 = 0
C(2,2) = 1; % u2 = 0
% Apoyo simple en 18
C(5,52) = 1; % u52 = 0
C(6,53) = 1; % u53 = 0
% Carretón en nodo 6 (impide giro y desplazamiento vertical pero está inclinado)
C(4,51) = 1; % u51 = 0
C(3,49) = -sin(pi/4); 
C(3,50) = cos(pi/4);  % -u49 + u50 = 0
% Desplazamiento impuesto en el nodo 4
C(7,32) = 1; % u32 = u (hacia abajo)

% MPCs

% Funciones lagrangianas paracada elemento ==> Diapositiva 56

x_1 = 5*L/16; x_2 = 6*L/16; x = L/3;
chi_mpc = (2/(x_2-x_1))*(x-(x_2+x_1)/2); % Diapositiva 47

n1 = (1-chi_mpc)/2; n2 = (1+chi_mpc)/2;

N_axil = [n1 n2];
% Flecha y giro desacoplados (independientes)
N_flecha = [n1 0 n2 0];
N_giro = [0 n1 0 n2];

C(8,16) = n1; C(8,19) = n2; % Desplazamiento x en la viga ppal
C(8,64) = -1; % Desplazamiento x en la riostra

C(9,17) = n1; C(9,20) = n2; % Desplazamiento y en la viga ppal
C(9,65) = -1; % Desplazamiento y en la riostra

C(10,18) = n1; C(10,21) = n2; % Giro en la viga ppal
C(10,66) = -1; % Giro en la riostra

C_quest = C;

%% Multiplicadores de Lagrange (No sirve porque la K sale singular)

KmL = [K C'; C zeros(p)];
r = zeros(p,1);
r(7) = u;           % Desplazamiento vertical en el nodo 4/11 en gdl 11/32, su posición es la línea 7
 
% FmL = [F; r]; % Vector de fuerzas ampliadas
% UmL = KmL\FmL;      % Aquí me avisa de que KmL es casi singular y los resultados pueden ser malos,
%                     % se vio en clase que esto podía pasar
% 
% u_mL = UmL([1:GDL]);
% lambda = UmL([GDL+1:end]); % Desde 27 hasta el final
% 
% f_R = -lambda;
% 
% f_uload = -lambda(7);

%% Método de la penalización

k_p = eye(p);   % Matriz diagonal de pxp
kpen = max(max(KmL))*1E6;       % Coeficiente de penalización
k = kpen*k_p;

Kpen = K + C'*k*C;
Fpen = F + C'*k*r;

Upen = Kpen\Fpen;

%% Cálculo de tensiones (Grader 4.1)

% Sobre la viga (creo un vector de desplazamientos que contenga los 18 gdl de la viga)

% u_ppal = [upen([1 2 3 4 5 6]) upen([4 5 6 7 8 9]) upen([7 8 9 10 11 12]) upen([10 11 12 13 14 15]) upen([13 14 15 16 17 18])];
% 
% chi = 0;
% B_f = 2*[0 -(1/2) 0 (1/2)];                         % El *2 es el iJe, se ralla el grader y por eso lo pongo a mano
% B_c = [-(1/2)*2 -(1-chi)/2 (1/2)*2 -(1+chi)/2];
% 
% for i=1:5 % 5 columnas de la matriz de desplazamientos de la viga ppal
% sigma_cort(i) = G*(B_c*u_ppal([2 3 5 6],i));
% sigma_fl(i) =  (-a1/2)*E*B_f*u_ppal([2 3 5 6],i);
% end
% 
% tau_max_v1 = max(abs(sigma_cort))/1e6;
% sigmaf_max_v1 = max(abs(sigma_fl))/1e6;
% 
% % Sobre la riostra (creo un vector de desplazamientos que contenga los 6 gdl de la riostra)
% 
% u_r1 = u_mL([19 20 21 22 23 24]);
% u_r2 = u_mL([22 23 24 25 26 27]);
% u_r = [u_r1 u_r2];
% 
% for i=1:2
% ur(:,i) = R*u_r(:,i);
% B = [-1 1]/Le_r;
% sigma_axil(i) = E*B*ur([1,4],i);
% end
% 
% sigma_axil_v2 = max(sigma_axil)/1e6;

%% Gráfica (Grader 4.2)

uv1 = Upen([1:GDL_v1]);
U_x1 = Upen([1:3:GDL_v1]); 
U_y1 = Upen([2:3:GDL_v1]);

% Deformada elástica

figure(1)

escala1 = 0.1*L/(max(abs(U_y1))); 
coord1x = coord_vl(:,1); coord1y = coord_vl(:,2);
grafu_1 = plot(coord1x,coord1y,'b*-',coord1x+escala1*U_x1,coord1y+escala1*U_y1,'r*--','LineWidth',2);
hold on

% Riostra
uv2 = Upen([GDL_v1+1:GDL_v1+GDL_r]);
u_x2 = Upen([GDL_v1+1:3:GDL_v1+GDL_r]);
u_y2 = Upen([GDL_v1+2:3:GDL_v1+GDL_r]); 

coordx2 = coord_r(:,1); coordy2 = coord_r(:,2);
grafu_2 = plot(coordx2,coordy2,'b*-',coordx2+escala1*u_x2,coordy2+escala1*u_y2,'gs--','LineWidth',2);

% Cálculo de desplazamientos

KmL = [K C'; C zeros(p)] ; 
F = zeros(size(K,1),1);
r = zeros(p,1); 
r(7) = u; 
FmL = [F; r]; 
UmL = KmL\FmL;

u_mL = UmL([1:GDL]); 
lambda = UmL([GDL+1:end]); 

f_R = -lambda;

f_uload = -lambda(7);
Le_pl = Le;
Je_pl = Je; 
iJe_pl = 6.4; 

x_i(:,:) = [Le/2:Le:L-Le/2]; 

u_ppal = [u_mL([1 2 3 4 5 6])       u_mL([4 5 6 7 8 9])       u_mL([7 8 9 10 11 12])    u_mL([10 11 12 13 14 15])...
          u_mL([13 14 15 16 17 18]) u_mL([16 17 18 19 20 21]) u_mL([19 20 21 22 23 24]) u_mL([22 23 24 25 26 27])...
          u_mL([25 26 27 28 29 30]) u_mL([28 29 30 31 32 33]) u_mL([31 32 33 34 35 36]) u_mL([34 35 36 37 38 39])...
          u_mL([37 38 39 40 41 42]) u_mL([40 41 42 43 44 45]) u_mL([43 44 45 46 47 48]) u_mL([46 47 48 49 50 51])];
      
for i=1:ne_v1
    % Momento flector
    B_f = iJe_pl*[0 -1/2 0 1/2];
    Mze_i(i,:) = E*I1*B_f*u_ppal([2 3 5 6],i); % Momento flector de los elementos en los puntos de integración en una matriz: en filas los elementos, en columnas los puntos de integracion
    % Fuerza cortante
    B_c = iJe_pl*[-1/2 0 1/2 0]-(1/2)*[0 1 0 1];
    AG = ((1/(Aef_v1*G))+((Le_pl^2)/(12*E*I1)))^-1;
    Vy_i(i,:) =  AG*B_c*u_ppal([2 3 5 6],i); % Fuerza cortante de flexion de los elementos en los puntos de integracion
end

figure(2)
grafmomment_int = plot(x_i(:,:),Mze_i(:,:)','bs-','LineWidth',1);
title('Momento flector Nm','Fontsize',10);
xlabel('posiciones nodales X','Fontsize',10) % Etiqueta el eje horizontal
ylabel('Nm','Fontsize',10) % Etiqueta el eje vertical
figure(4)
grafsforce_int = plot(x_i(:,:),Vy_i(:,:)','bs-','LineWidth',1);
title('Fuerza cortante N','Fontsize',10);
xlabel('posiciones nodales X','Fontsize',10) % Etiqueta el eje horizontal
ylabel('N','Fontsize',10) % Etiqueta el eje vertical


%% Funciones

% Matriz de Timoshenko, integración selectiva ===> Diapositiva 65 (misma notación)

function K_beam = Kbeam(L,I,A,E,G)
AG = ((1/(A*G))+((L^2)/(12*E*I)))^-1;
% Matriz de flexión
Kfe_i1 = [0 0 0 0; 0 E*I/L 0 -E*I/L; 0 0 0 0; 0 -E*I/L 0 E*I/L];
% Matriz de cortaruda
Kse_i1 = [AG/L  AG/2 -AG/L AG/2;  AG/2 AG*L/4 -AG/2 AG*L/4; -AG/L -AG/2 AG/L -AG/2; AG/2 AG*L/4 -AG/2 AG*L/4];
% Matriz de rigidez a flexuión y cortadura
K_beam = Kfe_i1 + Kse_i1;
end
