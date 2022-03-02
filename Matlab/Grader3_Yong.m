clear, close all
%% Datos de entrada
%
E = 70E9; nu = 0.3; L = 1; a0 = L/25; t = a0/25; q0 = 2200; P = -250; %%%%% OJO QUE YA LLEVA EL SIGNO
%

ne = 4; % número de elementos

nn = 5; % número de nodos

gdln = 3; %grados de libertad por nodo

coord_n = zeros(nn,2); % Matriz de coordenadas nodales: en filas nodos, en columna1 coordenada x, en columna2 coordenada y
                       % Numero de filas = numero de nodos, 2 columnas
coord_n(:,1) = 0:L/ne:L;
%
conectividad_e = zeros(ne,2);% Matriz de conectividad de elementos a través de los nodos
conectividad_e(:,1) = 1:1:ne;
conectividad_e(:,2) = 2:1:nn;



% solucion analitica



% Datos preliminares
%
ne = 4; % Numero de elementos (enunciado)
nn = ne+1; % Numero de nodos
gdln = 3; % GDL por nodo (3 porque es una viga)

% Datos de integracion numerica (copiado del archivo de clase: Clase_26_10_21)
%
chi_ip1 = 0; w_ip1 = 2; % cuadratura n=1
chi_ip2 = [-1 1]*1/sqrt(3); w_ip2 = [1 1]; % cuadratura n=2
chi_ip3 = [-sqrt(3/5) 0 sqrt(3/5)]; w_ip3 = [5/9 8/9 5/9]; % cuadratura n=3
%
%% Bucle de integracion y montaje de elementos
%
K = zeros(nn*gdln); % Declara la matriz de rigidez global
F = zeros(nn*gdln,1); % Declara el vector de fuerzas nodales en c.globales
%
for e = 1:ne; % Almacen de los grados de libertad de cada elemento   
    index = conectividad_e(e,:); % Indice de posicion de los nodos en global
    x1 = coord_n(index(1),1); x2 = coord_n(index(2),1);
    Le = x2-x1; Je = Le/2; iJe = 1/Je; % Jacobiano de la transformacion de cada elemento
    gdle = [index(1)*gdln-2 index(1)*gdln-1 index(1)*gdln...
            index(2)*gdln-2 index(2)*gdln-1 index(2)*gdln]; %indince de cada elemento
    gdlea =  gdle([1,4]); % Grados axiles
    gdlef =  gdle([2,3,5,6]); % Grados flectores
    %
    % Matriz de rigidez axil
    %
    Kea = zeros(2); 
    for i = 1:1 % Integracion numerica de Kea
        % Copiado y modificado del archivo de clase (Clase_19_10_21)
        chi = chi_ip1(i); w = w_ip1(i);
        Be = iJe*[-1/2 1/2]; % Para 2 nodos
        Ne = [(1-chi)/2 (1+chi)/2]; % Para 2 nodos
        %
        Ae_ = [A_(x1,a0,L,t); A_(x2,a0,L,t)]; %%%%%%%%%%%%%%%%%%%%%%%%%%% 
        Ae = Ne*Ae_;
        %
        Kea = Kea + Be'*E*Ae*Be*Je*w;
    end
    K(gdlea,gdlea)= K(gdlea,gdlea)+Kea; % Actualizacion de la rigidez axil en la matriz global
    %
    % Matriz de rigidez a flexion
    %
    Kef = zeros(4);
    %
    for i = 1:2 % Integracion numerica de Kef
        chi = chi_ip2(i); w = w_ip2(i);
        Ne = [(1-chi)/2 (1+chi)/2]; % Para 2 nodos
        % Copiado del archivo de clase (Clase_26_10_21)
        k1 = (chi*3/2); k2 = (chi*3-1)/2 ; k3 = -(chi*3/2); k4 = (chi*3+1)/2; % Derivada de __
        Bef = iJe^2*[k1 Je*k2 k3 Je*k4]; % Matriz cinematica de flexion
        %
        Ie_ = [Iy_(x1,a0,L,t);Iy_(x2,a0,L,t)]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%???????????????
        Ie = Ne*Ie_;
        %
        Kef = Kef + Bef'*E*Ie*Bef*Je*w;
        %
    end
    %
    K(gdlef,gdlef)=K(gdlef,gdlef)+Kef; % Actualizacion de la rigidez de flexion en la matriz global
    %
%     % Fuerzas axiles

    fea = zeros(2,1);
    %
    for i = 1:2 % Integracion numerica de fea
        chi = chi_ip2(i); w = w_ip2(i);
        Ne = [(1-chi)/2 (1+chi)/2]; % Para 2 nodos
        fea_ = [0;0];
        fea = fea+Ne'*Ne*fea_*Je*w;
    end
    %Da 0
    
%     %
    F(gdlea)= F(gdlea)+fea; % Actualizacion de la rigidez axil en la matriz global
%     %
    % Fuerzas de flexion
%     %
    fef = zeros(4,1);
%     %
    for i = 1:3 % Integracion numerica de fef
        chi = chi_ip3(i); w = w_ip3(i);
        NL = [(1-chi)/2 (1+chi)/2]; % Para 2 nodos;
        nh1 = -0.75*chi+0.25*chi^3+0.5; nh2 = -0.25*chi-0.25*chi^2+0.25*chi^3+0.25;
        nh3 = 0.75*chi-0.25*chi^3+0.5; nh4 = -0.25*chi+0.25*chi^2+0.25*chi^3-0.25;
        Neh = [nh1 Je*nh2 nh3 Je*nh4]; % Funciones hermiticas del elemento (de forma)
        
        fef_= [f_(L*3/4,f0,L);f_(L*7/8,f0,L)]; 
        
        fef = fef+Neh'*NL*fef_*Je*w;
    end
    F(gdlef) = F(gdlef)+fef; % Actualizacion de la rigidez de flexion en la matriz global
end
%
K_quest = K;
%
% %% Condiciones de contorno
% %
% u = zeros(nn*gdln,1);
% gdl_L = [4:nn*gdln]; gdl_R = [1 2 3]; % Grados libres (gdl_L) y restringidos (gdl_R)
% %
% K_LL = K(gdl_L,gdl_L);
% K_LR = K(gdl_L,gdl_R);
% %
% f_L = F(gdl_L); % Cargas nodos libres
% u_L = K_LL\f_L;
% %
% u(gdl_L) = u_L; % Vector de desplazamientos total
% %
% f_R = K_LR'*u_L; % Reacciones
% %
%% Funciones

function [xi_ip,w_ip] =Puntos_Pesos_Integracion(n)
if n==3
    xi_ip =[-sqrt(3/5), 0, sqrt(3/5)]; w_ip = [5/9, 8/9, 5/9];
elseif n==2    
    xi_ip = 1/sqrt(3)*[-1 1]; w_ip = [1 1];
end
end
 % if n==1    
%     xi_ip = 0; w_ip = 2;

% elseif n==3    
%     xi_ip =[-sqrt(3/5), 0, sqrt(3/5)]; w_ip = [5/9, 8/9, 5/9];
% elseif n==4    
%     xi_ip = [-sqrt((3+2*sqrt(6/5))/7),-sqrt((3-2*sqrt(6/5))/7), ...
%         sqrt((3-2*sqrt(6/5))/7), sqrt((3+2*sqrt(6/5))/7)];    
%     w_ip = [(18-sqrt(30))/36 , (18+sqrt(30))/36, ...
%         (18+sqrt(30))/36,(18-sqrt(30))/36];
%  elseif n==5   
%     xi_ip = [- 1/3 *sqrt(5+2*sqrt(10/7)),- 1/3 *sqrt(5-2*sqrt(10/7)),0, ...
%      1/3 *sqrt(5-2*sqrt(10/7)),  1/3 *sqrt(5+2*sqrt(10/7))];    
%     w_ip = [(322-13*sqrt(70))/900, (322+13*sqrt(70))/900,128/225,...        
%      (322+13*sqrt(70))/900,(322-13*sqrt(70))/900];
% else
%     disp('Error de input, la cuadratura maxima es 5')
% end

%
function Iy = Iy_(x,a0,L,t) % Inercia de la seccion
    a = a0*((-2/(3*L))*x+1); % Canto variable
    b = 3*a; % Ancho variable
    Iy = (1/12)*b*a^3-(1/12)*(b-2*t)*(a-2*t)^3; % Inercia de la seccion (respecto la linea media)
end
%

%Calculo de la Area ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function A = A_(x,a0,L,t) % Area de la seccion
    a = a0*((-2/(3*L))*x+1); % Canto variable
    b = 3*a; % Ancho variable
    A = a*b-(a-2*t)*(b-2*t); % Area de la seccion
%     A = 2*t*(a+b-2*t); % Area de la seccion    
end

% function A_9=A_(x,a0,L)
% A_9 = (a0*exp(-x/(3*L)))^2;%Area) a x a
% end

function Ne = N_3_nod(chi) %Funcion de forma  
Ne = [chi*(chi-1)/2 1-chi^2 chi*(chi+1)/2];
end

function Ne = N_2_nod(chi) %Funcion de forma  
Ne = [(1-chi)/2 (1+chi)/2];
end

%La Fuerza Calculada++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
function f=f_(x,f0,L)% la fuerza f(x) Calculado 
f = f0*(x/L)^4;
end   
%La Fuerza Calculada++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
function Be = B_3_nod(chi,Je)    
Be = (1/Je)*[chi-1/2 -2*chi chi+1/2];    
end    
    
%
function q = q_(x,q0,L) % Carga distribuida
    q = q0*(1-(x/L)^2);
end