% Grader 1 MEF

L = 1; q0 = 12000; E = 1000; I = 10; A = 100;
EI = E*I;
nele = 8; GDL = 24; gdln = 3;

%% Matriz de rigidez global

% Elemento 1

Le1 = L*sqrt(2)/2;
Ke1 = [E*A/Le1 0 0 -E*A/Le1 0 0; 0 12*EI/Le1^3 6*EI/Le1^2 0 -12*EI/Le1^3 6*EI/Le1^2; 0 6*EI/Le1^2 4*EI/Le1 0 -6*EI/Le1^2 2*EI/Le1 ; -E*A/Le1 0  0 E*A/Le1 0 0; 0 -12*EI/Le1^3 -6*EI/Le1^2 0 12*EI/Le1^3 -6*EI/Le1^2;  0 6*EI/Le1^2 2*EI/Le1 0 -6*EI/Le1^2 4*EI/Le1];    
alfa1 = pi/4;
R1 = [cos(alfa1) sin(alfa1) 0; -sin(alfa1) cos(alfa1) 0; 0 0 1];
zerosR = zeros(3);
R1_ = [R1 zerosR; zerosR R1];
Ke1_ = R1_'*Ke1*R1_;

% Elemento 2 (SE PUEDE OMITIR POR SER COMO 1)

Le2 = L*sqrt(2)/2;
Ke2 = [E*A/Le2 0 0 -E*A/Le2 0 0; 0 12*EI/Le2^3 6*EI/Le2^2 0 -12*EI/Le2^3 6*EI/Le2^2; 0 6*EI/Le2^2 4*EI/Le2 0 -6*EI/Le2^2 2*EI/Le2 ; -E*A/Le2 0  0 E*A/Le2 0 0; 0 -12*EI/Le2^3 -6*EI/Le2^2 0 12*EI/Le2^3 -6*EI/Le2^2;  0 6*EI/Le2^2 2*EI/Le2 0 -6*EI/Le2^2 4*EI/Le2];   
alfa2 = pi/4;
R2 = [cos(alfa2) sin(alfa2) 0; -sin(alfa2) cos(alfa2) 0; 0 0 1];
zerosR = zeros(3);
R2_ = [R2 zerosR; zerosR R2];
Ke2_ = R2_'*Ke2*R2_;

% Elemento 3

Le3 = L/2;
Ke3 = [E*A/Le3 0 0 -E*A/Le3 0 0; 0 12*EI/Le3^3 6*EI/Le3^2 0 -12*EI/Le3^3 6*EI/Le3^2; 0 6*EI/Le3^2 4*EI/Le3 0 -6*EI/Le3^2 2*EI/Le3 ; -E*A/Le3 0  0 E*A/Le3 0 0; 0 -12*EI/Le3^3 -6*EI/Le3^2 0 12*EI/Le3^3 -6*EI/Le3^2;  0 6*EI/Le3^2 2*EI/Le3 0 -6*EI/Le3^2 4*EI/Le3];    
alfa3 = 0;
R3 = [cos(alfa3) sin(alfa3) 0; -sin(alfa3) cos(alfa3) 0; 0 0 1];
zerosR = zeros(3);
R3_ = [R3 zerosR; zerosR R3];
Ke3_ = R3_'*Ke3*R3_;

% Elemento 4 (SE PUEDE OMITIR POR SER COMO 3)

Le4 = L/2;
Ke4 = [E*A/Le4 0 0 -E*A/Le4 0 0; 0 12*EI/Le4^3 6*EI/Le4^2 0 -12*EI/Le4^3 6*EI/Le4^2; 0 6*EI/Le4^2 4*EI/Le4 0 -6*EI/Le4^2 2*EI/Le4 ; -E*A/Le4 0  0 E*A/Le4 0 0; 0 -12*EI/Le4^3 -6*EI/Le4^2 0 12*EI/Le4^3 -6*EI/Le4^2;  0 6*EI/Le4^2 2*EI/Le4 0 -6*EI/Le4^2 4*EI/Le4];   
alfa4 = 0;
R4 = [cos(alfa4) sin(alfa4) 0; -sin(alfa4) cos(alfa4) 0; 0 0 1];
zerosR = zeros(3);
R4_ = [R4 zerosR; zerosR R4];
Ke4_ = R4_'*Ke4*R4_;

% Elemento 5

Le5 = L/2;
Ke5 = [E*A/Le5 0 0 -E*A/Le5 0 0; 0 12*EI/Le5^3 6*EI/Le5^2 0 -12*EI/Le5^3 6*EI/Le5^2; 0 6*EI/Le5^2 4*EI/Le5 0 -6*EI/Le5^2 2*EI/Le5 ; -E*A/Le5 0  0 E*A/Le5 0 0; 0 -12*EI/Le5^3 -6*EI/Le5^2 0 12*EI/Le5^3 -6*EI/Le5^2;  0 6*EI/Le5^2 2*EI/Le5 0 -6*EI/Le5^2 4*EI/Le5];    
alfa5 = -pi/2;
R5 = [cos(alfa5) sin(alfa5) 0; -sin(alfa5) cos(alfa5) 0; 0 0 1];
zerosR = zeros(3);
R5_ = [R5 zerosR; zerosR R5];
Ke5_ = R5_'*Ke5*R5_;

% Elemento 6 (SE PUEDE OMITIR POR SER COMO 5)

Le6 = L/2;
Ke6 = [E*A/Le6 0 0 -E*A/Le6 0 0; 0 12*EI/Le6^3 6*EI/Le6^2 0 -12*EI/Le6^3 6*EI/Le6^2; 0 6*EI/Le6^2 4*EI/Le6 0 -6*EI/Le6^2 2*EI/Le6 ; -E*A/Le6 0  0 E*A/Le6 0 0; 0 -12*EI/Le6^3 -6*EI/Le6^2 0 12*EI/Le6^3 -6*EI/Le6^2;  0 6*EI/Le6^2 2*EI/Le6 0 -6*EI/Le6^2 4*EI/Le6];    
alfa6 = -pi/2;
R6 = [cos(alfa6) sin(alfa6) 0; -sin(alfa6) cos(alfa6) 0; 0 0 1];
zerosR = zeros(3);
R6_ = [R6 zerosR; zerosR R6];
Ke6_ = R6_'*Ke6*R6_;

% Elemento 7

Le7 = L*sqrt(2)/2;
Ke7 = [E*A/Le7 0 0 -E*A/Le7 0 0; 0 12*EI/Le7^3 6*EI/Le7^2 0 -12*EI/Le7^3 6*EI/Le7^2; 0 6*EI/Le7^2 4*EI/Le7 0 -6*EI/Le7^2 2*EI/Le7 ; -E*A/Le7 0  0 E*A/Le7 0 0; 0 -12*EI/Le7^3 -6*EI/Le7^2 0 12*EI/Le7^3 -6*EI/Le7^2;  0 6*EI/Le7^2 2*EI/Le7 0 -6*EI/Le7^2 4*EI/Le7];    
alfa7 = 3*pi/4;
R7 = [cos(alfa7) sin(alfa7) 0; -sin(alfa7) cos(alfa7) 0; 0 0 1];
zerosR = zeros(3);
R7_ = [R7 zerosR; zerosR R7];
Ke7_ = R7_'*Ke7*R7_;

% Elemento 8 (SE PUEDE OMITIR POR SER COMO 7)

Le8 = L*sqrt(2)/2;
Ke8 = [E*A/Le8 0 0 -E*A/Le8 0 0; 0 12*EI/Le8^3 6*EI/Le8^2 0 -12*EI/Le8^3 6*EI/Le8^2; 0 6*EI/Le8^2 4*EI/Le8 0 -6*EI/Le8^2 2*EI/Le8 ; -E*A/Le8 0  0 E*A/Le8 0 0; 0 -12*EI/Le8^3 -6*EI/Le8^2 0 12*EI/Le8^3 -6*EI/Le8^2;  0 6*EI/Le8^2 2*EI/Le8 0 -6*EI/Le8^2 4*EI/Le8];  
alfa8 = 3*pi/4;
R8 = [cos(alfa8) sin(alfa8) 0; -sin(alfa8) cos(alfa8) 0; 0 0 1];
zerosR = zeros(3);
R8_ = [R8 zerosR; zerosR R8];
Ke8_ = R8_'*Ke8*R8_;

% Ensamblaje de los elementos == Sumatorio de fuerzas internas

K_quest = zeros(24);

% Elemento 1 
    
gdle1 = [1 2 3 4 5 6];
K_quest(gdle1, gdle1) = K_quest(gdle1, gdle1) + Ke1_;

% Elemento 2 

gdle2 = [4 5 6 7 8 9];
K_quest(gdle2, gdle2) = K_quest(gdle2, gdle2) + Ke2_;

% Elemento 3 

gdle3 = [7 8 9 10 11 12]; 
K_quest(gdle3, gdle3) = K_quest(gdle3, gdle3) + Ke3_;

% Elemento 4

gdle4 = [10 11 12 13 14 15];
K_quest(gdle4, gdle4) = K_quest(gdle4, gdle4) + Ke4_;
    
% Elemento 5 

gdle5 = [13 14 15 16 17 18];
K_quest(gdle5, gdle5) = K_quest(gdle5, gdle5) + Ke5_;

% Elemento 6

gdle6 = [16 17 18 19 20 21];
K_quest(gdle6, gdle6) = K_quest(gdle6, gdle6) + Ke6_;
    
% Elemento 7

gdle7 = [19 20 21 22 23 24];
K_quest(gdle7, gdle7) = K_quest(gdle7, gdle7) + Ke7_;
    
% Elemento 8

gdle8 = [22 23 24 7 8 9];                       % OJO
K_quest(gdle8, gdle8) = K_quest(gdle8, gdle8) + Ke8_;

K_quest = round(K_quest);

%% Fuerzas nodales equivalentes (Resolviendo la elástica en los elementos)

feq_quest = zeros(24,1);

f1x = 0; f1y = - q0*Le1/15; f1z = - q0*Le1^2/60;                                % Fuerzas nodales en el nodo 1 en coordenadas locales
f2x = 0; f2y = -4*q0*Le1/15 - q0*Le2/2; f2z = -q0*Le1^2/12 + q0*Le2^2/30;
f3x = 0; f3y = -q0*Le2/2; f3z = q0*Le2^2/12;

fn1 = [f1x f1y f1z]; fn1_ = R1'*fn1';                                           % Fuerzas nodales en coordenadas locales
fn2 = [f2x f2y f2z]; fn2_ = R1'*fn2';
fn3 = [f3x f3y f3z]; fn3_ = R1'*fn3';

feq_quest(1:9) = [fn1_ fn2_ fn3_];                                              % Vector completo de fuerzas nodales equivalentes


%% Condensación estática (impongo CC)

gdl_L = [4:18 22:24];
KLL = K_quest(gdl_L, gdl_L);                  % Matriz K reducida con CC en gdl_L
gdl_R = [1:3 19:21];
KRR = K_quest(gdl_R, gdl_R);                  % Matriz K reducida con CC en gdl_R
KRL = K_quest(gdl_R, gdl_L);
KLR = KRL';

fL = feq_quest(gdl_L,1);                      % Vector de fuerzas nodales en los nodos libres
uL = KLL\fL;                                  % Desplazamientos en los nodos centrales
gdl_uquest = [1:3 7:9 13:18];
u_quest = uL(gdl_uquest,1);

fR = KRL*uL;                                  % Vector de fuerzas de reacción en los apoyos
fr_quest = fR(1:3) - feq_quest(1:3);          % Fuerzas de reacción en el apoyo izquierdo (nodo 1) RESTO LAS NODALES PORQUE NO SON DE REACCIÓN


%% Verificación de resultados

u = zeros(24,1);
u(gdl_L) = uL;
f = K_quest*u;                                % Sus componentes deben ser las de fL y fR

%% Apartado 5

u7 = u(19:24);                                % Desplazamientos globales en el elemento 7
u7_local = R7_*u7;                            % Desplazamientos locales en el elemento 7

fele_qest =  Ke7*u7_local;                    % Fuerzas internas en el apoyo 7 == Las calculo con la matriz local y los desplazamientos locales

fele_ = f (19:24);
fele_quest = R7_*fele_;

%% Representación gráfica de cosas chulas
% 
% clf
% escala = 0.25/max(max(abs(ux),abs(uy)));
% figure(1)
% %whitebg(figure(1), [1 1 1])
% 
% for e=1:nele
%     indice = conect_ele(e,:)
%     gdl = gdle(e,:);
%     ue = u(gdl);   % Desplazamientos de cada elemento
%     n1 = coord_nod(indice(1),:); n2 = coord_nod(indice(2),:); conod_ele = [n1;n2];      % Coordenadas nodales de cada elemento
%     plot(conod_ele(:,1),conod_ele(:,2),'b*-', conod_ele(:,1)+escala*ue(1:3:6),conod_ele(:,2)+escala*ue(2:3:6),'r*--')
%     hold on
% end
% grid on
% axis equal
% title('Deformada elástica en autoescala 25%', 'Fontsize', 14)
% xlabel('Posiciones nodales X', 'Fontsize',12)
% ylabel('Posiciones nodales Y', 'Fontize', 12)

% figure(2)
% % whitebg(figure(1), [1 1 1])
% verK = zeros(GDL);
% for s=1:size(K_quest,1)
%     for t=1:size(K_quest,1)
%         if abs(K_quest(s,t))>1e-12
%             versK(s,t)=1;
%         else versK(s,t)<1e-12;
%         end
%     end
% end
% imagesc(verK);                  % Crea un plot de colores con los valores de K
% colormap(flipud(gray));         % Cambia el mapa de colores a gris
%                                 % los mayores negro y los menores gris
% title('Ocupación de K','Fontsize',14)
% xlabel('gdln','Fontsize',12)
% ylabel('gdln','Fontsize',12)
%            