%% GRADER 1: Repaso de Cálculo Matricial

% Rafael Rivero de Nicolás
% MEF

clc
clear

%% Enunciado, datos:

L_ = 1; q0 = 12000; E = 1000; I = 10; A = 100;

% CÁLCULO DE PARÁMETROS
L(1) = sqrt(2)*L_/2;    % Elemento 1 (Nodos 1-2)
alpha(1) = +pi/4;

L(2) = sqrt(2)*L_/2;    % Elemento 2 (Nodos 2-3)
alpha(2) = +pi/4;

L(3) = L_/2;    % Elemento 3 (Nodos 3-4)
alpha(3) = 0;

L(4) = L_/2;    % Elemento 4 (Nodos 4-5)
alpha(4) = 0;

L(5) = L_/2;    % Elemento 5 (Nodos 5-6)
alpha(5) = -pi/2;

L(6) = L_/2;    % Elemento 6 (Nodos 6-7)
alpha(6) = -pi/2;

L(7) = sqrt(2)*L_/2;    % Elemento 7 (Nodos 7-8)
alpha(7) = +3*pi/4;

L(8) = sqrt(2)*L_/2;    % Elemento 8 (Nodos 8-3)
alpha(8) = +3*pi/4;

%% Definicicón de las matrices k_ en ejes globales

for i = 1:8
    k(:,:,i) = [E*A/L(i)    0           0         -E*A/L(i)      0                0;       ...      %Matriz de rigidez en locales (6x6)
                0   12*E*I/L(i)^3    6*E*I/L(i)^2     0      -12*E*I/L(i)^3   6*E*I/L(i)^2;...
                0   6*E*I/L(i)^2     4*E*I/L(i)       0      -6*E*I/L(i)^2    2*E*I/L(i);  ...
                -E*A/L(i)   0          0           E*A/L(i)      0                0;       ...
                0   -12*E*I/L(i)^3   -6*E*I/L(i)^2    0      12*E*I/L(i)^3   -6*E*I/L(i)^2;...
                0   6*E*I/L(i)^2     2*E*I/L(i)       0      -6*E*I/L(i)^2    4*E*I/L(i)];
            
            
    R(:,:,i) = [cos(alpha(i)) sin(alpha(i)) 0; -sin(alpha(i)) cos(alpha(i)) 0; 0 0 1]; %(3x3)

    R_(:,:,i) = [R(:,:,i) zeros(3); zeros(3) R(:,:,i)]; %Matriz de rotación ampliada
    
    k_(:,:,i) = (R_(:,:,i))'*k(:,:,i)*R_(:,:,i); % Matriz de rigidez en ejes globales
end   

%% Ensamblaje de la Matriz de rigidez de la estructura en ejes globales

K = zeros(24); %Inicialización de la matriz de rigidez de la estructura

gdl(1,:) = [1 2 3 4 5 6]; %El elemento 1 afecta a los grados de libertad 1, 2, 3, 4, 5, 6.
gdl(2,:) = [4 5 6 7 8 9];
gdl(3,:) = [7 8 9 10 11 12];
gdl(4,:) = [10 11 12 13 14 15];
gdl(5,:) = [13 14 15 16 17 18];
gdl(6,:) = [16 17 18 19 20 21];
gdl(7,:) = [19 20 21 22 23 24];
gdl(8,:) = [22 23 24 7 8 9];

for i = 1:8 
    
    K(gdl(i,:),gdl(i,:)) = K(gdl(i,:),gdl(i,:)) + k_(:,:,i);
    
end 

K = round(K);

% % % % Lo que hace el bucle:
% % % % K(gdl1, gdl1) = K(gdl1, gdl1) + ke1_ %Ensamblaje del elemento 1
% % % % K(gdl2, gdl2) = K(gdl2, gdl2) + ke2_ %Ensamblaje del elemento 2
% % % % K(gdl3, gdl3) = K(gdl3, gdl3) + ke3_ %Ensamblaje del elemento 3

%% Fuerzas nodales

% Nodo 1
Fn(1,:) =  R(:,:,1)'*( [0; -q0*L(1)/15; (-q0*L(1)^2)/60] );

% Nodo 2 
Fn(2,:) =  R(:,:,1)'*( [0; -4*q0*L(1)/15; q0*L(1)^2/30] + [0; -q0*L(2)/2; -q0*L(2)^2/12] );

% Nodo 3
Fn(3,:) = R(:,:,2)'*( [0; -q0*L(2)/2; q0*L(2)^2/12] );

%Nodo 1
feq_quest(1) = Fn(1,1);
feq_quest(2) = Fn(1,2);
feq_quest(3) = Fn(1,3);
%Nodo 2
feq_quest(4) = Fn(2,1);
feq_quest(5) = Fn(2,2);
feq_quest(6) = Fn(2,3);
%Nodo 3
feq_quest(7) = Fn(3,1);
feq_quest(8) = Fn(3,2);
feq_quest(9) = Fn(3,3);

for j = 10:24
    feq_quest(j) = 0;
end

%F_L(1:gdl_L) FL = DATO 

% %Nodo 2
FL(1:3,1) = feq_quest(4:6);
% % FL(1) = ;%gdl 4
% % FL(2) = ;%gdl 5
% % FL(3) = ;%gdl 6
% %Nodo 3
FL(4:6,1) = feq_quest(7:9);
% % FL(4) = ;%gdl 7
% % FL(5) = ;%gdl 8
% % FL(6) = ;%gdl 9

% % Nodos 4, 5, 6 & 8
for j = 7:18
    FL(j,1) = 0;
end

%% Condensación Estática 

%Grados de libertad libres
gdl_L = [4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 22 23 24]; 

%Grados de libertad restringidos
gdl_R = [1 2 3 19 20 21]; %(Apoyos)

%Matriz de rigidez reducida por las condiciones de contorno
KRR = K(gdl_R,gdl_R); %(6x6)
KLL = K(gdl_L,gdl_L); %(18x18)
KRL = K(gdl_R,gdl_L); %(6x6)
KLR = KRL';

% % % % % [------------------]
% % % % % |   KRR  |  KRL    |
% % % % % | -------|-------- |
% % % % % |   KLR  |  KLL    |
% % % % % [------------------]

uL = KLL\FL; % inv(KLL)
FR = KRL*uL;

u = [0 0 0 uL(1:15)' 0 0 0 uL(16:18)'];
f = [FR(1:3)' FL(1:15,1)' FR(4:6)' FL(16:18)'];

%% Verificación 

u_p = zeros(18,1);
u_p(gdl_L) = uL;
f_p = K*u_p;
Verif = f'-f_p;% Debería ser prácticamente cero

%% Desplazamientos de los nodos 2, 4, 6, 8

u_quest = zeros(12,1);
u_quest = [u(4:6), u(10:12), u(16:18), u(22:24)]';

%% Fuerzas de reacción en el Nodo 1

fr_quest = [FR(1:3)] - [Fn(1,1) Fn(1,2) Fn(1,3)]';

%% Fuerza Interna Elemento 7

u_7 = R_(:,:,7)*[u(19:24)]';
fele_quest = k(:,:,7)*u_7;

%% SOLUCIONES DEL GRADER 

% % % % % % vector de fuerzas nodales equivalentes en c. globales
% % % % % feq_quest = feq_quest'
% % % % % q_quest = feq_quest;
% % % % % % Matriz de rigidez global
% % % % % 
% % % % % K_quest = K
% % % % % 
% % % % % % desplazamiento del nodo central
% % % % % 
% % % % % u_quest
% % % % % 
% % % % % % fuerza de reacción en apoyo izquierdo
% % % % % f_quest = fele_quest













