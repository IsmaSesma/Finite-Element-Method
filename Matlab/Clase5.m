%% Ejercicio MPS barra y viga

%% Datos
EI = 1; L = 1; EAv = 100; P = -1;

% Elemento viga 1
Le1 = L;
Kviga11 = [EAv/Le1 0 0; 0 12*EI/Le1^3 6*EI/Le1^2; 0 6*EI/Le1^2 4*EI/Le1];
Kviga12 = [-EAv/Le1 0 0; 0 -12*EI/Le1^3 6*EI/Le1^3; 0 -6*EI/Le1^2 2*EI/Le1];
Kviga22 = [EAv/Le1 0 0; 0 12*EI/Le1^3 -6*EI/Le1^2; 0 -6*EI/Le1^2 4*EI/Le1];
Kviga = [Kviga11 Kviga12 ; Kviga12' Kviga22];


%% Elemento barra (slave en este caso)
Le2 = L*sqrt(2)/2;
EAb = EAv/10;
Kbar = EAb/Le2*[1 0 -1 0;0 0 0 0;-1 0 1 0; 0 0 0 0];
alfa = pi/4;
Rb = [cos(alfa) sin(alfa); -sin(alfa) cos(alfa)];
Rb_ = [Rb zeros(2); zeros(2) Rb];
Kbarra = Rb_'*Kbar*Rb_;


%% Rigidez del modelo
zeros_ = zeros(length(Kviga),length(Kbarra));
K = [Kviga zeros_; zeros_' Kbarra];


%% Matriz C
r = zeros(7,1);
p = length(r);
C = zeros(length(r), size(K,1));
% SPCs
C(1,1) = 1;
C(2,2) = 1;
C(3,3) = 1;
C(4,7) = 1;
C(5,8) = 1;
% MPCs
chi_mpc = 0;
n1 = (1-chi_mpc)/2; n2 =(1+chi_mpc)/2;
C(6,1) = n1; C(6,4) = n2; C(6,9) = -1;
Je1 = Le1/2;

nh1 = -0.75*chi_mpc + 0.25*chi_mpc^3 + 0.5;
nh2 = Je1*(-0.25*chi_mpc - 0.25*chi_mpc^2 + 0.25*chi_mpc^3 + 0.25);
nh3 = 0.75*chi_mpc - 0.25*chi_mpc^3 + 0.5;
nh4 = Je1*(-0.25*chi_mpc + 0.25*chi_mpc^2 + 0.25*chi_mpc^3 - 0.25);

C(7,2) = nh1; C(7,3) = nh2; C(7,5) = nh3; C(7,6) = nh4; C(7,10) = -1;

C;

KmL = [K C'; C zeros(p)]; % matriz aumentada
F = zeros(size(K,1),1);
F(5) = P;
FmL = [F; r]; % vector fuerzas ampliadas
UmL = KmL\FmL;

u_mL = UmL([1:10])
lambda = UmL([11:end]) % los dos últimos términos son las fuerzas internas que se transmiten entre los elementos
f_R = -lambda

% para comprobar si está bien calculas la inversa de la KmL


%% Penalización
k_p = eye(p); % matriz diagonal de pxp
kpen = max(max(KmL))*10^3;% coeficiente de penalización
k = kpen*k_p;
Kpen = K + C'*k*C;
Fpen = F + C'*k*r;
Upen = Kpen\Fpen;
error = Upen-u_mL;







