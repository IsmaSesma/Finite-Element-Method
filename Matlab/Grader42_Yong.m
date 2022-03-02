%% Datos
%
E = 70E9; nu = 0.3;
G = E/(2*(1+nu)); 
L = 5; 
%
a1 = L/40; 
b1 = a1/2;
t1 = a1/10; 
kap_1 = 2.5; 
%
r2 = a1/5; 
t2 = r2/5; 
kap_2 = 2; 
%
A1 = 2*t1*b1+t1*(a1-2*t1); 
A2 = pi*r2^2-pi*(r2-t2)^2; 
%
nn_e = 2; 
ne_pl = 16;
nn_pl = ne_pl+1; 
ne_r = 4; 
nn_r = ne_r+1; 
%
u = -L/100; 
%
Iv1 = (1/12)*t1*(a1-2*t1)^3+2*((1/12)*b1*t1^3+b1*t1*((a1-t1)/2)^2);
Iv2 = (pi/4)*(r2^4-(r2-t2)^4);
%
%integracion numerica
%
chi_ip1 = 0; w_ip1 = 2; 
chi_ip2 = [-1 1]*1/sqrt(3); w_ip2 = [1 1]; 
chi_ip3 = [-sqrt(3/5) 0 sqrt(3/5)]; w_ip3 = [5/9 8/9 5/9];
%
% Coordenada de los nodos
%
coordn_pl = zeros(nn_pl,2); 
coordn_pl(:,1) = 0:L/ne_pl:L; 
coordn_r = zeros(nn_r,2); 
coordn_r(:,1) = 0:(L/3)/ne_r:L/3;
coordn_r(:,2) = -L/2:(L/2)/ne_r:0;
%
%% Bucles de integracion y montaje de los elementos
%
% Grados de libertad
%
gdln_pl = 3;
GDL_pl = nn_pl*gdln_pl;
gdln_r = 3; 
GDL_r = nn_r*gdln_r; 
GDL = GDL_pl + GDL_r; 
%
% Inicializacion de las variables
%
K = zeros(GDL); 
F = zeros(GDL,1);

%viga ppal 
%
conectividad_e = zeros(ne_pl,2); 
conectividad_e(:,1) = 1:1:ne_pl;
conectividad_e(:,2) = 2:1:nn_pl;
%
for e = 1:ne_pl %
    coord_n = coordn_pl;
    gdln = gdln_pl;
    %
    index = conectividad_e(e,:); % Indice de posicion de los nodos en global
    x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
    Le = x2-x1; Je = Le/2; iJe = 1/Je; % Jacobiano de la transformacion de cada elemento
    gdle = [index(1)*gdln-2 index(1)*gdln-1 index(1)*gdln...
            index(2)*gdln-2 index(2)*gdln-1 index(2)*gdln];
    gdlea =  gdle([1,4]); % Grados axiles
    gdlef =  gdle([2,3,5,6]); % Grados flectores
   
    % Matriz de rigidez axil
   
    Kea_pl = zeros(2); % Inicializacion
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
        Kea_pl = Kea_pl + Be'*E*Ae*Be*Je*w;
        %
    end
    %
    K(gdlea,gdlea) = K(gdlea,gdlea) + Kea_pl; % Actualizacion de la rigidez axil en la matriz global
    %
end
%
A_pl = A1/kap_1; % 
Kefc_pl = Kbeam(Le,Iv1,A_pl,E,G);
%
% Ensamblaje
%
Kefc = zeros(GDL); % Inicializacion (establecer el tamaño)
%
for i = 2:3:47;
    Kefc([i i+1 i+3 (i+3)+1],[i i+1 i+3 (i+3)+1]) = Kefc([i i+1 i+3 (i+3)+1],[i i+1 i+3 (i+3)+1]) + Kefc_pl;
end
%
K = K + Kefc; 

%  viga enrriostrada 

Le_r = sqrt((L/3)^2+(L/2)^2)/ne_r; % Longitud de cada elemento de la viga enrriostrada (Longitud total de la viga enrriostrada/numero de elementos de la viga enrriostrada)
Je_r = Le_r/2; % Jacobiano
iJe_r = Je^-1; % Inversa del jacobiano
%
A_r = A2/kap_2; % Area efectiva de la viga enrriostrada
%
Kefc_r = Kbeam(Le_r,Iv2,A_r,E,G);
%
% Matriz de rigidez axil
%
Kea_r = zeros(nn_e*gdln_r); % Inicializacion
%
for i=1:1
    %
    Je = Je_r;
    chi = chi_ip1(i); w = w_ip1(i);
    Be = B_(chi,Je);
    Ne = N_(chi);
    %
    Ae_ = [A2 A2];
    Ae = Ne*[A2 A2]'; % Traspuesta para que las dimensiones sean correctas
    %
    Kea_r([1 4],[1 4]) = Kea_r([1 4],[1 4]) + Be'*E*Ae*Be*Je*w; % [1,4] porque son los gdl axiles
    %
end
%
Ke_r_ = zeros(nn_e*gdln_r); % Matriz rigidez de un elemento de la viga enrriostrada
Ke_r_([2 3 5 6],[2 3 5 6]) = Ke_r_([2 3 5 6],[2 3 5 6]) + Kefc_r;
Ke_r = Kea_r + Ke_r_;

% Cambio de base

alpha = atan((L/2)/(L/3));
K_e_r = R_(alpha)'*Ke_r*R_(alpha); 

% Ensamblaje de la matriz de rigidez 

K_f_r = zeros(nn_r*gdln_r);
%
for i = 1:3:10
    K_f_r([i i+1 i+2 i+3 i+4 i+5],[i i+1 i+2 i+3 i+4 i+5]) = K_f_r([i i+1 i+2 i+3 i+4 i+5],[i i+1 i+2 i+3 i+4 i+5]) + K_e_r;
end

% Ensamblaje de la matriz de rigidez global 

K([GDL_pl+1:GDL],[GDL_pl+1:GDL]) = K([GDL_pl+1:GDL],[GDL_pl+1:GDL]) + K_f_r; % Para que a lo que ya estaba en la matriz K (procedente de la viga ppal) se le sume de lo de la riostra en esas componentes

K_quest = K;

%% Calculo de la matriz de restricciones
%
p = 10; 
C = zeros(p,size(K,1)); 
% SPCs


C(1,1) = 1;
C(2,2) = 1; 

C(3,49) = -sin(pi/4); C(3,50) = cos(pi/4); 
C(4,51) = 1; 

C(5,52) = 1;
C(6,53) = 1; 
% Desplazamiento impuesto entre los nodos 10 y 11
x_1 = 9*L/16; x_2 = 10*L/16; x = 0.6*L+0.125;
chi = (2/(x_2-x_1))*(x-(x_2+x_1)/2);
N1 = (1/2)*(1-chi); N2 = (1/2)*(1+chi);
C(7,29) = N1; C(7,32) = N2;
%
% MPCs

x_1 = 5*L/16; x_2 = 6*L/16; x = L/3;
chi = (2/(x_2-x_1))*(x-(x_2+x_1)/2); 
%
N1 = (1/2)*(1-chi); N2 = (1/2)*(1+chi);
%
N_a = [N1 N2];
%
N_f = [N1 0 N2 0];
N_g = [0 N1 0 N2];
%
C(8,16) = N1; C(8,19) = N2; 
C(8,64) = -1; 
%
C(9,17) = N1; C(9,20) = N2;
C(9,65) = -1; 
%
C(10,18) = N1; C(10,21) = N2;
C(10,66) = -1; 

C_quest = C;
%
%% Penalizacion 
%
KmL = [K C'; C zeros(p)] ; 
r = zeros(p,1); 
r(7) = u; 

k_p = eye(p);
kpen = 1E6*max(max(KmL));
k = kpen*k_p;
%
Kpen = K + C'*k*C;
%
Fpen = F+C'*k*r;
%
Upen = Kpen\Fpen;

%% Gráfica 

uv1 = Upen([1:GDL_pl]);
U_x1 = Upen([1:3:GDL_pl]); 
U_y1 = Upen([2:3:GDL_pl]);

% deformada elástica

figure(1)

escala1 = 0.1*L/(max(abs(U_y1))); 
coord1x = coordn_pl(:,1); coord1y = coordn_pl(:,2);
grafu_1 = plot(coord1x,coord1y,'b*-',coord1x+escala1*U_x1,coord1y+escala1*U_y1,'r*--','LineWidth',2);
hold on

% la riostra
uv2 = Upen([GDL_pl+1:GDL_pl+GDL_r]);
u_x2 = Upen([GDL_pl+1:3:GDL_pl+GDL_r]);
u_y2 = Upen([GDL_pl+2:3:GDL_pl+GDL_r]); 
%
coordx2 = coordn_r(:,1); coordy2 = coordn_r(:,2);
grafu_2 = plot(coordx2,coordy2,'b*-',coordx2+escala1*u_x2,coordy2+escala1*u_y2,'gs--','LineWidth',2);
%

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
iJe_pl = iJe; 
%
x_i(:,:) = [Le/2:Le:L-Le/2]; 
%
u_ppal = [u_mL([1 2 3 4 5 6])       u_mL([4 5 6 7 8 9])       u_mL([7 8 9 10 11 12])    u_mL([10 11 12 13 14 15])...
          u_mL([13 14 15 16 17 18]) u_mL([16 17 18 19 20 21]) u_mL([19 20 21 22 23 24]) u_mL([22 23 24 25 26 27])...
          u_mL([25 26 27 28 29 30]) u_mL([28 29 30 31 32 33]) u_mL([31 32 33 34 35 36]) u_mL([34 35 36 37 38 39])...
          u_mL([37 38 39 40 41 42]) u_mL([40 41 42 43 44 45]) u_mL([43 44 45 46 47 48]) u_mL([46 47 48 49 50 51])];
for i=1:ne_pl
    % Momento flector
    B_f = iJe_pl*[0 -1/2 0 1/2];
    Mze_i(i,:) = E*Iv1*B_f*u_ppal([2 3 5 6],i); % Momento flector de los elementos en los puntos de integración en una matriz: en filas los elementos, en columnas los puntos de integracion
    % Fuerza cortante
    B_c = iJe_pl*[-1/2 0 1/2 0]-(1/2)*[0 1 0 1];
    AG = ((1/(A_pl*G))+((Le_pl^2)/(12*E*Iv1)))^-1;
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
%
%% Funciones

function N = N_(chi) 
    N = [(1-chi)/2 (1+chi)/2];
end
%
function B = B_(chi,Je) %
    
    iJe = Je^-1;
    B = iJe*[-1/2 1/2];
end
%
function d2_N_h = d2_N_h_(chi,Je) 
    iJe = Je^-1;
    d2_nh1 = 3*chi/2; 
    d2_nh2 = (-1+3*chi)/2; 
    d2_nh3 = -3*chi/2; 
    d2_nh4 = (1+3*chi)/2;
    d2_N_h = iJe^2 * [d2_nh1 Je*d2_nh2 d2_nh3 Je*d2_nh4]; 
end
%
function K_beam = Kbeam(L,I,A,E,G) 
    AG = ((1/(A*G))+((L^2)/(12*E*I)))^-1;
   
    Kfe_i1 = [0 0 0 0;
              0 E*I/L 0 -E*I/L;
              0 0 0 0;
              0 -E*I/L 0 E*I/L];
    
    Kse_i1 = [AG/L AG/2 -AG/L AG/2; 
              AG/2 AG*L/4 -AG/2 AG*L/4;
             -AG/L -AG/2 AG/L -AG/2;
              AG/2 AG*L/4 -AG/2 AG*L/4];
    
    K_beam = Kfe_i1 + Kse_i1;
end

function R = R_(alpha)
    R_ = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1];
    R = [R_ zeros(3); zeros(3) R_];
end