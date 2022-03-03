% **************************************** 
%        Lesson 0 FEM Matrix Calculus
% ****************************************

clc;clear;close all;

E = 100; A = 1; L = 1; I = 100;
EI = E*I;

%% Element 1

Le1 = L*sqrt(2);
Ke1 = [E*A/Le1 0 0 -E*A/Le1 0 0,...
    0 12*EI/Le1^3 6*EI/Le1^2 0 -12*EI/Le1^3 6*EI/Le1^2,...
    0 6*EI/Le1^2 4*EI/Le1 0 -6*EI/Le1^2 2*EI/Le1,...
    -E*A/Le1 0  0 E*A/Le1 0 0,...                                                       % Stiffness matrix
    0 -12*EI/Le1^3 -6*EI/Le1^2 0 12*EI/Le1^3 -6*EI/Le1^2,...
    0 6*EI/Le1^2 2*EI/Le1 0 -6*EI/Le1^2 4*EI/Le1];     

alfa1 = -pi/4;
R1 =  [cos(alfa1) sin(alfa1) 0; -sin(alfa1) cos(alfa1) 0; 0 0 1];                       % Rotation matrix for bar and beam 
zerosR = zeros(3);
R1_ = [R1 zerosR; zerosR R1]; % Base change matrix 
Ke1_ = R1_'*Ke1*R1_;


%% Element 2

Le2 = L;
Ke2 = E*A/Le2*[1 0 -1  0; 0 0 0 0;  -1 0 1 0; 0 0 0 0]; % Needs to be updated to take into account the beam's DOF
alfa2 = 0;
R2 =  [cos(alfa2) sin(alfa2) 0; -sin(alfa2) cos(alfa2) 0; 0 0 1];
zerosR = zeros(3);
R2_ = [R2 zerosR; zerosR R2];
Ke2_ = R2_'*Ke2*R2_; 

%% Element 3

Le3 = L*sqrt(1.5^2 + 1^2);
Ke3 = E*A/Le3*[1 0 -1 0; 0 0 0 0;  -1 0 1 0; 0 0 0 0];  % Needs to be updated to take into account the beam's DOF
alfa3 = atan(1.5);
R3 =  [cos(alfa3) sin(alfa3) 0; -sin(alfa3) cos(alfa3) 0; 0 0 1];
zerosR = zeros(3);
R3_ = [R3 zerosR; zerosR R3];
Ke3_ = R3_'*Ke3*R3_; 

%% Element assembly (internal forces sumation)

K = zeros(12);

% Element 1

gdle1 = [1 2  3 10 11 12];
K(gdle1, gdle1) = K(gdle1, gdle1) + Ke1_;   % Ensamblaje del elemento 1

% Element 2

gdle2 = [4 5 6 10 11 12];
K(gdle2, gdle2) = K(gdle2, gdle2) + Ke2_;   % Ensamblaje del elemento 2

% Element 3

gdle3 = [7 8 9 10 11 12];
K(gdle3, gdle3) = K(gdle3, gdle3) + Ke3_;   % Ensamblaje del elemento 3


%% Boundary conditions (aplication)

gdl_L = [3 6 9 10 11 12];
KLL = K(gdl_L, gdl_L);              % Stiffness matrix with boundary conditions in free DOF (gdl_L)
gdl_R = [1:9];
KRR = K(gdl_R, gdl_R);              % Stiffness matrix with boundary conditions in restricted DOF (gdl_R)
KRL = K(gdl_R, gdl_L);
KLR = KRL';

F = 100;
alfa = pi/4;
fL = F*[0 0 0 cos(alfa) -sin(alfa) 0];     % External forces vector in free DOF
uL = KLL\fL;                               % uL = KLL^(-1)*fL is slower and wont be used

fR = KRL*uL;                               % Reaction forces in restricted nodes 
%% Verificationn

u = zeros(12,1);
u(gdl_L)= uL;    % Displacement vector in all DOF 
f = K*u;         % Force vector (its component must match fR and fL)
