clear, close all
%% Datos
%General
E = 70E9; nu = 0.3;
G = E/(2*(1+nu)); % Modulo de cizalladura de un material isotropo
L = 5; % L viga principal

%Grado de Liberdad definiendo
nn_e = 2; % Numero nodos por elemento
nn_pl = 6; % Numero nodos viga principal ( hay nodos: 1-6)
ne_pl = 5; % Numero elementos viga de Timoshenko en la viga principal
nn_r = 3; % Numero nodos viga riostra (hay nodos: 7-9)
ne_r = 2; % Numero elementos viga de Timoshenko en la viga riostra 
%
%Viga Principal
a1 = L/40; % Canto de la seccion recta
b1 = a1/2; % Ancho de la seccion recta
t1 = a1/10; % Espesor del perfil
kap_1 = 2.5; % Coeficiente de cortante para el area de la viga
%Riostra
r2 = a1/5; % Radio 2 de la riostra
t2 = r2/5; % Espesor tubular de la riostra
kap_2 = 2; % Coeficiente de cortante para el area de la riostra
%Areas
A2 = pi*r2^2-pi*(r2-t2)^2; % Area de la viga arriostrada
A1 = 2*t1*b1+t1*(a1-2*t1); % Area de la viga ppal

%Condicion restrigido

u = -L/100; % Condicion de desplazmiento vertical impuesto en el nodo 4 (el signo porque va hacia abajo)

%Momento de inercia
Iv1 = (1/12)*t1*(a1-2*t1)^3+2*((1/12)*b1*t1^3+b1*t1*((a1-t1)/2)^2); % Momento de inercia de la viga ppal
Iv2 = (pi/4)*(r2^4-(r2-t2)^4); % Momento de inercia de la viga enriostrada

% integracion numerica 

chi_ip1 = 0; w_ip1 = 2; % Cuadratura n=1
chi_ip2 = [-1 1]*1/sqrt(3); w_ip2 = [1 1]; % Cuadratura n=2
chi_ip3 = [-sqrt(3/5) 0 sqrt(3/5)]; w_ip3 = [5/9 8/9 5/9]; % Cuadratura n=3
%
% Coordenada de los nodos
%
coordn_pl = zeros(nn_pl,2); % Coordenadas de los nodos de la viga ppal
% En las columnas define: 1 es x,  2 es y
% hay 6 nodos de la viga ppal
coordn_pl(:,1) = 0:L/ne_pl:L; % Va De 0 hasta L (0:L), con salto de L/nn_pl
% esta inicialmente 0 en y 
coordn_r = zeros(nn_r,2); % Coordenadas de los nodos de la viga enrriostrada
%
coordn_r(:,1) = 0:(L/3)/ne_r:L/3; % en x salta de L/3/ne_r
coordn_r(:,2) = -L/2:(L/2)/ne_r:0;%en y mide -L/2
%
%%  integracion y montaje de elementos 
%
% Grados de libertad
%
gdln_pl = 3; % GDL por nodo de la viga pl
GDL_pl = nn_pl*gdln_pl; % GDL totales de la viga pl
gdln_r = 3; % Grados de libertad por nodo de la vida enrriostrada
GDL_riostra = nn_r*gdln_r; % Grados de libertad totales de la viga enrriostrada
GDL = GDL_pl + GDL_riostra; % Grados totales de la estructura
%
% Inicializacion de las variables
%
K = zeros(GDL); % Inicializacion de la matriz de rigidez global
F = zeros(GDL,1); % Inicializacion del vector de fuerzas nodales en coordenadas globales
%
% Calculo de la viga ppal (viga de Timoshenko)
%
conectividad_e = zeros(ne_pl,2); % Matriz de conectividad de elementos a través de los nodos
conectividad_e(:,1) = 1:1:ne_pl;
conectividad_e(:,2) = 2:1:nn_pl;
%
for e = 1:ne_pl % Almacen de los grados de libertad de cada elemento
%
% Copia y pega del resto de codigos (practicamente)
%
coord_n = coordn_pl;
gdln = gdln_pl;
% 
index = conectividad_e(e,:); % Indice de posicion de los nodos en global
x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
Le = x2-x1; Je = Le/2; iJe = 1/Je; % Jacobiano de la transformacion de cada elemento
gdle = [index(1)*gdln-2 index(1)*gdln-1 index(1)*gdln...
index(2)*gdln-2 index(2)*gdln-1 index(2)*gdln];
gdlea = gdle([1,4]); % Grados axiles
gdlef = gdle([2,3,5,6]); % Grados flectores
%
% Matriz de rigidez axil
%
Kea_ppal = zeros(2); % Inicializacion
%
for i=1:1
%
chi = chi_ip1(i); w = w_ip1(i);
Be = B_(chi,Je);
Ne = N_(chi);
%
Ae_ = [A1;A1]; % Area gorrito
Ae = Ne*Ae_; % Area virgulilla
%
Kea_ppal = Kea_ppal + Be'*E*Ae*Be*Je*w;
%
end
%
K(gdlea,gdlea) = K(gdlea,gdlea) + Kea_ppal; % Actualizacion de la rigidez axil en la matriz global
%
end
%
A_ppal = A1/kap_1; % Area efectiva de la viga ppal
%
Kefc_ppal = Kbeam(Le,Iv1,A_ppal,E,G);
%
% Ensamblaje
%
Kefc = zeros(GDL); % Inicializacion (establecer el tamaño)
%
Kefc([2 3 5 6],[2 3 5 6]) = Kefc([2 3 5 6],[2 3 5 6]) + Kefc_ppal;
Kefc([5 6 8 9],[5 6 8 9]) = Kefc([5 6 8 9],[5 6 8 9]) + Kefc_ppal;
Kefc([8 9 11 12],[8 9 11 12]) = Kefc([8 9 11 12],[8 9 11 12]) + Kefc_ppal;
Kefc([11 12 14 15],[11 12 14 15]) = Kefc([11 12 14 15],[11 12 14 15]) + Kefc_ppal;
Kefc([14 15 17 18],[14 15 17 18]) = Kefc([14 15 17 18],[14 15 17 18]) + Kefc_ppal;
%
K = K + Kefc; % A lo que habia ya en el K (la rigidez axil aportada por la viga ppal), se la suma lo de la matriz de flexion+cortadura

%
% Calculo de la viga enrriostrada (viga Timoshenko)
%
Le_riostra = sqrt((L/3)^2+(L/2)^2)/ne_r; % Longitud de cada elemento de la viga enrriostrada (Longitud total de la viga enrriostrada/numero de elementos de la viga enrriostrada)
Je_riostra = Le_riostra/2; % Jacobiano
iJe_riostra = Je^-1; % Inversa del jacobiano
%
A_riostra = A2/kap_2; % Area efectiva de la viga enrriostrada
%
Kefc_riostra = Kbeam(Le_riostra,Iv2,A_riostra,E,G);
%
% Matriz de rigidez axil
%
Kea_riostra = zeros(nn_e*gdln_r); % Inicializacion
%
for i=1:1
%
Je = Je_riostra;
chi = chi_ip1(i); w = w_ip1(i);
Be = B_(chi,Je);
Ne = N_(chi);
%
Ae_ = [A2 A2];
Ae = Ne*[A2 A2]'; % Traspuesta para que las dimensiones sean correctas
%
Kea_riostra([1 4],[1 4]) = Kea_riostra([1 4],[1 4]) + Be'*E*Ae*Be*Je*w; % [1,4] porque son los gdl axiles
%
end
%
Ke_riostra_ = zeros(nn_e*gdln_r); % Matriz rigidez de un elemento de la viga enrriostrada
Ke_riostra_([2 3 5 6],[2 3 5 6]) = Ke_riostra_([2 3 5 6],[2 3 5 6]) + Kefc_riostra;
Ke_riostra = Kea_riostra + Ke_riostra_;
%
% Cambio de base de la matriz de rigidez de la viga enrriostrada
%
alpha = atan((L/2)/(L/3));
K_e_riostra = R_(alpha)'*Ke_riostra*R_(alpha); % Matriz rigidez de un elemento de la viga enrriostrada ya rotada
%
% Ensamblaje de la matriz de rigidez de la viga enrriostrada (las aportaciones de rigidez de cada elemento a la viga enrriostrada)
%
K_final_riostra = zeros(nn_r*gdln_r);
K_final_riostra([1 2 3 4 5 6],[1 2 3 4 5 6]) = K_final_riostra([1 2 3 4 5 6],[1 2 3 4 5 6]) + K_e_riostra; % La rigidez que aporta el elemento 1 de la viga enrriostrada
K_final_riostra([4 5 6 7 8 9],[4 5 6 7 8 9]) = K_final_riostra([4 5 6 7 8 9],[4 5 6 7 8 9]) + K_e_riostra; % La rigidez que aporta el elemento 2 de la viga enrriostrada
%
% Ensamblaje de la matriz de rigidez global (lo que aporta la viga enrriostrada en conjunto a la matriz de rigidez global de la estructura)
%
K([GDL_pl+1:GDL],[GDL_pl+1:GDL]) = K([GDL_pl+1:GDL],[GDL_pl+1:GDL]) + K_final_riostra; % Para que a lo que ya estaba en la matriz K (procedente de la viga ppal) se le sume de lo de la riostra en esas componentes
%
K_quest = K;
%
%% Calculo de la matriz de restricciones
%
% p = 9; % Numero de restricciones
p = 10; % Numero de restricciones (9 de apoyos + 1 desplazamiento impuesto)
%
C = zeros(p,size(K,1)); % Declaracion de la matriz de coeficientes de restriccion
%
% SPCs
%
% Apoyo simple en nodo 1
C(1,1) = 1; % u1 = 0
C(2,2) = 1; % u2 = 0
% Apoyo simple en nodo 7
C(5,19) = 1; % u19 = 0
C(6,20) = 1; % u20 = 0
% Carriton en nodo 6 (impide giro y desplazamiento en su direccion vertical, en este caso esta inclinado)
C(4,18) = 1; % u18 = 0
C(3,16) = -sin(pi/4); C(3,17) = cos(pi/4);% -u16 + u17 = 0
% Desplazamiento impuesto en el nodo 4
C(7,11) = 1; % u11 = u
%
% MPCs
%
% Funciones lagrangianas para elemento viga de Timoshenko en 1D (Diapositiva 56)

x_1 = L/5; x_2 = 2*L/5; x = L/3;
chi = (2/(x_2-x_1))*(x-(x_2+x_1)/2); % Diapositiva 47

N1 = (1/2)*(1-chi); N2 = (1/2)*(1+chi);

N_axil = [N1 N2];
% Flecha y giro desacoplados
N_flecha = [N1 0 N2 0];
N_giro = [0 N1 0 N2];

C(8,4) = N1; C(8,7) = N2; % Desplazamiento x en la viga ppal
C(8,25) = -1; % Desplazamiento x en la riostra

C(9,5) = N1; C(9,8) = N2; % Desplazamiento y en la viga ppal
C(9,26) = -1; % Desplazamiento y en la riostra

C(10,6) = N1; C(10,9) = N2; % Giro en la viga ppal
C(10,27) = -1; % Giro en la riostra

C_quest = C;

%% Metodo de los multiplicadores de Lagrange
%
KmL = [K C'; C zeros(p)] ; % Matriz ampliada
%
F = zeros(size(K,1),1);
r = zeros(p,1); % Restricciones
r(7) = u; % Desplazamiento vertical en el nodo 4 (u11): La ecuacion esta en la linea 7
%
FmL = [F; r]; % Vector de fuerzas ampliadas
UmL = KmL\FmL;
%
% Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND = 1.170825e-16. 
% > In Untitled5 (line 229)
%
u_mL = UmL([1:GDL]); % Desde 1 hasta 27
lambda = UmL([GDL+1:end]); % Desde 27 hasta el final
%
f_R = -lambda;
%
f_uload = -lambda(7); % 
%
% %% Penalizacion (mirarlo, porque no me coincide Upen con u_mL)
% %
% k_p = eye(p);
% kpen = 1E6*max(max(KmL)); % Coeficiente de penalizacion
% % max(KmL): maximo de cada fila
% % El maximo de la diagonal principal es tbn el maximo de la matriz de rigidez, sino nos da, es que lo hemos hecho mal
% %
% k = kpen*k_p;
% %
% Kpen = K + C'*k*C;
% %
% Fpen = F+C'*k*r;
% %
% Upen = Kpen\Fpen;
% %
%% Tension axil sobre la riostra
%
u_e1 = u_mL([19 20 21 22 23 24]);
u_e2 = u_mL([22 23 24 25 26 27]);
u_riostra = [u_e1 u_e2];
%
for i=1:2
u_riostra_local(:,i) = R_(alpha)*u_riostra(:,i);
B = [-1 1]/Le_riostra;
sigma_axil(i) = E*B*u_riostra_local([1,4],i);
end
%
sigma_axil_v2 = max(sigma_axil)/1e6;
%
%% Tensiones de la viga ppal
%
u_ppal = [u_mL([1 2 3 4 5 6]) u_mL([4 5 6 7 8 9]) u_mL([7 8 9 10 11 12]) u_mL([10 11 12 13 14 15]) u_mL([13 14 15 16 17 18])];
%
chi = 0; % Punto de integración del elemento (Diapositiva 16: Integracion abierta -> Regla del pto medio: en el pto medio, chi es 0)
%
B_f = iJe*[0 -(1/2) 0 (1/2)];
B_c = [-(1/2)*iJe -(1-chi)/2 (1/2)*iJe -(1+chi)/2];
%
for i=1:5 % 5 columnas de la matriz de desplazamientos de la viga ppal
sigma_cortadura(i) = G*(B_c*u_ppal([2 3 5 6],i));
sigma_flexion(i) = (-a1/2)*E*B_f*u_ppal([2 3 5 6],i);
end
%
tau_max_v1 = max(abs(sigma_cortadura))/1e6;
sigmaf_max_v1 = max(abs(sigma_flexion))/1e6;


%% Funciones
%
function N = N_(chi) % Funcion de forma de 2 nodos
% N = [N1 N2];
N = [(1-chi)/2 (1+chi)/2];
end
%
function B = B_(chi,Je) %
% B = (Je)^-1*[N1' N2'];
iJe = Je^-1;
B = iJe*[-1/2 1/2];
end
%
function d2_N_h = d2_N_h_(chi,Je) % Segunda derivada funcion hermitica
iJe = Je^-1;
d2_nh1 = 3*chi/2; 
d2_nh2 = (-1+3*chi)/2; 
d2_nh3 = -3*chi/2; 
d2_nh4 = (1+3*chi)/2;
d2_N_h = iJe^2 * [d2_nh1 Je*d2_nh2 d2_nh3 Je*d2_nh4]; % Coincide con la matriz cinematica y de curvaturas (diapositiva 47)
end
%
function K_beam = Kbeam(L,I,A,E,G) % Diapositiva 65: Matriz viga de Timoshenko
% Matriz de rigidez a flexion (Diapositiva 62 - Igual para la integracion exacta y para integracion selectiva-reducida)
Kfe_i1 = [0 0 0 0;
0 E*I/L 0 -E*I/L;
0 0 0 0;
0 -E*I/L 0 E*I/L];
% Matriz de rigidez a cortadura (Diapositiva 62 - Integracion selectiva-reducida)
Kse_i1 = [A*G/L A*G/2 -A*G/L A*G/2; 
A*G/2 A*G*L/4 -A*G/2 A*G*L/4;
-A*G/L -A*G/2 A*G/L -A*G/2;
A*G/2 A*G*L/4 -A*G/2 A*G*L/4];
% Matriz de rigidez a flexion+cortadura (Diapositiva 65)
K_beam = Kfe_i1 + Kse_i1;
end
%
function R = R_(alpha)% Matriz de rotacion
R_ = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1];
R = [R_ zeros(3); zeros(3) R_];
end