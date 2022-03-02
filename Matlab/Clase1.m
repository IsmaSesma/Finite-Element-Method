%% Clase 0 de cálculo matricial MEF
E = 100; A = 1; L = 1; I = 100;
EI = E*I;

%% Elemento 1

Le1 = L*sqrt(2);
Ke1 = [E*A/Le1 0 0 -E*A/Le1 0 0; 0 12*EI/Le1^3 6*EI/Le1^2 0 -12*EI/Le1^3 6*EI/Le1^2; 0 6*EI/Le1^2 4*EI/Le1 0 -6*EI/Le1^2 2*EI/Le1 ; -E*A/Le1 0  0 E*A/Le1 0 0; 0 -12*EI/Le1^3 -6*EI/Le1^2 0 12*EI/Le1^3 -6*EI/Le1^2;  0 6*EI/Le1^2 2*EI/Le1 0 -6*EI/Le1^2 4*EI/Le1];    % Matriz de rigidez del elemento 1
alfa1 = -pi/4;
R1 =  [cos(alfa1) sin(alfa1) 0; -sin(alfa1) cos(alfa1) 0; 0 0 1];                       % Matriz de giro con barra y viga
zerosR = zeros(3);
R1_ = [R1 zerosR; zerosR R1]; % Matriz de cambio de base del elemento 1
Ke1_ = R1_'*Ke1*R1_;


%% Elemento 2

Le2 = L;
Ke2 = E*A/Le2*[1 0 -1 0; 0 0 0 0;  -1 0 1 0; 0 0 0 0];
alfa2 = 0;
R2 =  [cos(alfa2) sin(alfa2) 0; -sin(alfa2) cos(alfa2) 0; 0 0 1];
zerosR = zeros(3);
R2_ = [R2 zerosR; zerosR R2];
Ke2_ = R2_'*Ke2*R2_; % Matriz de rigidez del elemento 2 en coordenadas globales

%% Elemento 3

Le3 = L*sqrt(1.5^2 + 1^2);
Ke3 = E*A/Le3*[1 0 -1 0; 0 0 0 0;  -1 0 1 0; 0 0 0 0];
alfa3 = atan(1.5);
R3 =  [cos(alfa3) sin(alfa3) 0; -sin(alfa3) cos(alfa3) 0; 0 0 1];
zerosR = zeros(3);
R3_ = [R3 zerosR; zerosR R3];
Ke3_ = R3_'*Ke3*R3_; % Matriz de rigidez del elemento 3 en coordenadas globales

%% Ensamblaje de los elemntos (sumatorio de fuerzas internas)

K = zeros(12); % Declaración o inicialización de la matriz de rigidez

%% Elemento 1

gdle1 = [1 2  3 10 11 12];
K(gdle1, gdle1) = K(gdle1, gdle1) + Ke1_;   % Ensamblaje del elemento 1

%% Elemento 2

gdle2 = [4 5 6 10 11 12];
K(gdle2, gdle2) = K(gdle2, gdle2) + Ke2_;   % Ensamblaje del elemento 2

%% Elemento 3

gdle3 = [7 8 9 10 11 12];
K(gdle3, gdle3) = K(gdle3, gdle3) + Ke3_;   % Ensamblaje del elemento 3


%% Condiciones de contorno (aplicación)

gdl_L = [3 6 9 10 11 12];
KLL = K(gdl_L, gdl_L);              % Matriz de rigidez reducida por las condiciones de contorno en gdl_L
gdl_R = [1:9];
KRR = K(gdl_R, gdl_R);              % Matriz de rigidez reducida por las condiciones de contorno en gdl_R
KRL = K(gdl_R, gdl_L);
KLR = KRL';

F = 100;
alfa = pi/4;
fL = F*[0 0 0 cos(alfa) -sin(alfa) 0];     % Vector de fuerzas externas en grados libres
uL = KLL\fL;                               % uL = KLL^(-1)*fL es más lento y no lo usaremos

fR = KRL*uL;                               % Fuerzas de reacción en los nodos restringidos

%% Verificación

u = zeros(12,1);
u(gdl_L)= uL;    % Vector completo de desplazamientos en los gdl del problema
f = K*u;         % Vector completo de fuerzas (sus componentes tienen que coincidir con fR y fL)




