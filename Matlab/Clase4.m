% Problema de la viga hermítica de 2 nodos (grados flectores) // Clase 4

E = 1;  A = 100;  I = 1; L = 1; f0y = -1; Py = 10;

%% Datos Integración Numérica

%Cuadratura n=1
chi_ip1 = 0; %Coordenadas de los puntos de integración
w_ip1 = 2;   %Peso de integración
%Cuadratura n=2
chi_ip2 = [-1 1]*1/sqrt(3); %Coordenadas de los puntos de integración
w_ip2 = [1 1]; %Peso de integración
%Cuadratura n=3
chi_ip3 = [-sqrt(3/5) 0 sqrt(3/5)]; %Coordenadas de los puntos de integración
w_ip3 = [5/9 8/9 5/9]; %Peso de integración
%Cuadratura n=4
chi_ip4 = [-0.8611363116 -0.3399810436 0.3399810436 0.8611363116];
w_ip4 = [0.3478548451 0.6521451549 0.6521451549 0.3478548451];


ipk = 2; %Cuadratura de integración de K
if ipk==1
    chi_ipk = chi_ip1; w_ipk = w_ip1;
elseif ipk==2
    chi_ipk = chi_ip2; w_ipk = w_ip2;
elseif ipk==3
    chi_ipk = chi_ip3; w_ipk = w_ip3;
elseif ipk==4
    chi_ipk = chi_ip4; w_ipk = w_ip4;
end

ipf = 3; %Cuadratura de integración de Fuerzas
if ipf==1
    chi_ipf = chi_ip1; w_ipf = w_ip1;
elseif ipf==2
    chi_ipf = chi_ip2; w_ipf = w_ip2;
elseif ipf==3
    chi_ipf = chi_ip3; w_ipf = w_ip3;
elseif ipf==4
    chi_ipf = chi_ip4; w_ipf = w_ip4;
end


%% Matriz de rigidez

x2 = L; x1 = 0;
Le = (x2-x1);
Je = Le/2;
iJe = Je^(-1);

Kef = zeros(4);              % 4 DOF

for ip = 1:ipk                          % Si integro con orden 1 me queda una matriz diagonal que solo reproduce el giro y da lugar a un mecanismo (tendrá 0s en la diagonal)
                                        % Si integro con orden 3 el resultado va a ser el mismo
    chi = chi_ipk(ip);
    w = w_ipk(ip);
    k1 = (chi*3/2); k2 = (chi*3-1)/2; k3 = -(chi*3)/2; k4 = (chi*3+1)/2;
    Bef = iJe^2*[k1 Je*k2 k3 Je*k4];
    Kef = Kef + Bef'*E*I*Bef*Je*w;
end

%% Vector de fuerzas nodales equivalentes

Feqf = zeros(4,1);                            % Mismos DOF que desplazamientos
for ip = 1 : ipf
    chi = chi_ipf(ip); w = w_ipf(ip);

    nh1 = -0.75*chi+0.25*chi^3+0.5; nh2 = -0.25*chi-0.25*chi^2+0.25*chi^3+0.25;
    nh3 = 0.75*chi-0.25*chi^3+0.5;  nh4 = -0.25*chi+0.25*chi^2+0.25*chi^3-0.25;
    Nhe = [nh1 Je*nh2 nh3 Je*nh4];            % Funciones de forma hermíticas del elemento

    f0y_ = [0;f0y];
    NLe = [(1-chi)/2 (1+chi)/2];        % Funciones de formas lagrangianas

    f_chi = NLe*f0y_;
    Feqf = Feqf + Nhe'*f_chi*Je*w;      
    
end