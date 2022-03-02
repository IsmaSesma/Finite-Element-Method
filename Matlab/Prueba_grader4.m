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
nn_v1 = 6;          % Nodos de la primera viga
nn_r = 3;           % Nodos de la riostra
ne_v1 = 5;          % Número de elementos en la viga 1
ne_r = 2;           % Número de elementos de la riostra

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

conectividad_e = zeros(ne_v1,2); % Matriz de conectividad de elementos a través de los nodos
conectividad_e(:,1) = 1:1:ne_v1;
conectividad_e(:,2) = 2:1:nn_v1;

% Elemento viga 1

a1 = L/40; b1 = a1/2; t1 = a1/10;
kc1 = 2.5;      % Coeficiente de cortante en la viga 1
I1 = 1/12*t1*(a1 - 2*t1)^3 + 2*(1/12*b1*t1^3 + b1*t1*((a1 - t1)/2)^2);  % Momento de inercia de la viga 1
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

Kef([2 3 5 6],[2 3 5 6]) = Kef([2 3 5 6],[2 3 5 6]) + Kef_v1;
Kef([5 6 8 9],[5 6 8 9]) = Kef([5 6 8 9],[5 6 8 9]) + Kef_v1;
Kef([8 9 11 12],[8 9 11 12]) = Kef([8 9 11 12],[8 9 11 12]) + Kef_v1;
Kef([11 12 14 15],[11 12 14 15]) = Kef([11 12 14 15],[11 12 14 15]) + Kef_v1;
Kef([14 15 17 18],[14 15 17 18]) = Kef([14 15 17 18],[14 15 17 18]) + Kef_v1;
% Como Timoshenko no considera el axil lo añado
K = K + Kef;           % Añado a la rigidez axil la flexión/cortadura por Timoshenko

% Para la riostra (la considero también una viga)

Le_r = sqrt((L/3)^2+(L/2)^2)/ne_r;          % Longitud de cada elemento de la riostra
Je_r = Le_r/2;                              % Jacobiano de la transformación
iJe = Je_r^(-1);                            % Inversa del jacobiano

Aef_r = A2/kc2;                             % Área efectiva de la riostra

Kef_r = Kbeam(Le_r, I2, Aef_r, E, G);

% Matriz de rigidez axil de la riostra

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

Ke_r_f = zeros(GDL_r);
Ke_r_f([1 2 3 4 5 6],[1 2 3 4 5 6]) = Ke_r_f([1 2 3 4 5 6],[1 2 3 4 5 6]) + Ke_r_giro; % Elemento 1 de la riostra
Ke_r_f([4 5 6 7 8 9],[4 5 6 7 8 9]) = Ke_r_f([4 5 6 7 8 9],[4 5 6 7 8 9]) + Ke_r_giro; % Elemento 2 de la riostra

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
% Apoyo simple en 7
C(5,19) = 1; % u19 = 0
C(6,20) = 1; % u20 = 0
% Carretón en nodo 6 (impide giro y desplazamiento vertical pero está inclinado)
C(4,18) = 1; % u18 = 0
C(3,16) = -sin(pi/4); 
C(3,17) = cos(pi/4);  % -u16 + u17 = 0
% Desplazamiento impuesto en el nodo 4
C(7,11) = 1; % u11 = u (hacia abajo)

% MPCs

% Funciones lagrangianas paracada elemento ==> Diapositiva 56

x_1 = L/5; x_2 = 2*L/5; x = L/3;
chi_mpc = (2/(x_2-x_1))*(x-(x_2+x_1)/2); % Diapositiva 47

n1 = (1-chi_mpc)/2; n2 = (1+chi_mpc)/2;

N_axil = [n1 n2];
% Flecha y giro desacoplados (independientes)
N_flecha = [n1 0 n2 0];
N_giro = [0 n1 0 n2];

C(8,4) = n1; C(8,7) = n2; % Desplazamiento x en la viga ppal
C(8,25) = -1; % Desplazamiento x en la riostra

C(9,5) = n1; C(9,8) = n2; % Desplazamiento y en la viga ppal
C(9,26) = -1; % Desplazamiento y en la riostra

C(10,6) = n1; C(10,9) = n2; % Giro en la viga ppal
C(10,27) = -1; % Giro en la riostra

C_quest = C;

%% Multiplicadores de Lagrange

KmL = [K C'; C zeros(p)];
r = zeros(p,1);
r(7) = u;           % Desplzamiento vertical en el nodo 4 en gdl 11, su posición es la línea 7

FmL = [F; r]; % Vector de fuerzas ampliadas
UmL = KmL\FmL;      % Aquí me avisa de que KmL es casi singular y los resultados pueden ser malos,
                    % se vio en clase que esto podía pasar

u_mL = UmL([1:GDL]);
lambda = UmL([GDL+1:end]); % Desde 27 hasta el final

f_R = -lambda;

f_uload = -lambda(7);


%% Cálculo de tensiones

% Sobre la viga (creo un vector de desplazamientos que contenga los 18 gdl de la viga)

u_ppal = [u_mL([1 2 3 4 5 6]) u_mL([4 5 6 7 8 9]) u_mL([7 8 9 10 11 12]) u_mL([10 11 12 13 14 15]) u_mL([13 14 15 16 17 18])];

chi = 0;
B_f = iJe*[0 -(1/2) 0 (1/2)];
B_c = [-(1/2)*iJe -(1-chi)/2 (1/2)*iJe -(1+chi)/2];

for i=1:5 % 5 columnas de la matriz de desplazamientos de la viga ppal
sigma_cort(i) = G*(B_c*u_ppal([2 3 5 6],i));
sigma_fl(i) = (-a1/2)*E*B_f*u_ppal([2 3 5 6],i);
end

tau_max_v1 = max(abs(sigma_cort))/1e6;
sigmaf_max_v1 = max(abs(sigma_fl))/1e6;

% Sobre la riostra (creo un vector de desplazamientos que contenga los 6 gdl de la riostra)

u_r1 = u_mL([19 20 21 22 23 24]);
u_r2 = u_mL([22 23 24 25 26 27]);
u_r = [u_r1 u_r2];

for i=1:2
ur(:,i) = R*u_r(:,i);
B = [-1 1]/Le_r;
sigma_axil(i) = E*B*ur([1,4],i);
end

sigma_axil_v2 = max(sigma_axil)/1e6;


%% Funciones

% Matriz de Timoshenko, integración selectiva ===> Diapositiva 65 (misma notación)

function K_beam = Kbeam(L,I,A,E,G)
% Matriz de flexión
Kfe_i1 = [0 0 0 0; 0 E*I/L 0 -E*I/L; 0 0 0 0; 0 -E*I/L 0 E*I/L];
% Matriz de cortaruda
Kse_i1 = [A*G/L  A*G/2 -A*G/L A*G/2;  A*G/2 A*G*L/4 -A*G/2 A*G*L/4; -A*G/L -A*G/2 A*G/L -A*G/2; A*G/2 A*G*L/4 -A*G/2 A*G*L/4];
% Matriz de rigidez a flexuión y cortadura
K_beam = Kfe_i1 + Kse_i1;
end