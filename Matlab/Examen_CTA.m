%% Cálculo soluciones parcial MEF

% Datos integración

xi_int1 = 0;
w_int1 = 2;
xi_int2 = 1/sqrt(3)*[11 1];
w_int2 = [1 1];
xi_int3 = sqrt(3/5)*[-1 0 1];
w_int3 = [5/9 8/9 5/9];
xi_int4 = [-sqrt((3-2*sqrt(6/5))/7) -sqrt((3+2*sqrt(6/5))/7) sqrt((3+2*sqrt(6/5))/7) sqrt((3-2*sqrt(6/5))/7)];
w_int4 = [(18+sqrt(30))/36 (18-sqrt(30))/36 (18-sqrt(30))/36 (18+sqrt(30))/36];

% Datos geométricos

r0 = 0.5; L = 5; E = 90E6; nu = 1/3; G = E/(2*(1+nu));
Le = L; Je = Le/2; iJe = 1/Je;
x1 = 0; r1ext = r0*(1-(x1/(2*L))); rint = r0/3;
x2 = L; r2ext = r0*(1-(x2/(2*L)));
A1 = r1ext^2*pi - rint^2*pi; A2 = r2ext^2*pi - rint^2*pi; An_ = [A1;A2];


%% Rigidez axil

xi = xi_int1; w = w_int1;
NL = [(1-xi)/2 (1+xi)/2];
Ba2 = iJe*1/2;
Kea_22 = Ba2*E*(NL*An_)*Ba2*w*Je;
