% ***********************************************************************
%           FEM MODEL OF A PLATE USING RECTANGULAR ELEMENTS 
% ***********************************************************************

clc;clear;close all;

%% MASIC AND GEOMETRIC DATA

plate.E = 70E9; plate.nu = 0.3;
plate.a = 1; plate.b = 1;                           % Plate's dimensions
plate.t = 0.01;                                      % Plate's thickness
plate.m = 0.5;                                      % Plate's mass
plate.rho = plate.m/plate.a/plate.b/plate.t;        % Plate's density
plate.I = plate.t^3/12;                             % Plate's inertia
plate.G = plate.E/2/(1 + plate.nu);

%% NUMERIC INTEGRATION DATA (Gauss-Legendre)

% n = 1 quadrature
chi_ip.f = 0;   % Integration point coordinate
w_ip.f = 2;     % Integration weight
% n = 2 quadrature
chi_ip.s = [-1 1]*1/sqrt(3);   
w_ip.s = [1 1];    
% Integration order
iob = 2;     % Bending matrix  
ios = 1;     % Shear matrix
iof = 1;     % Force vector

%% INPUT DATA AND MESH CREATION 

ne_x = 10;                                              % Number of elements in X
ne_y = ne_x;                                            % Number of elements in Y
structure = CreateMesh(plate.a, plate.b,ne_x,ne_y);     % Mesh definition
ne = length(structure.mesh.elements.id);                % Number of elements
nne = 4;                                                % Number of nodes per element
nn = length(structure.mesh.nodes.id);                   % Number of nodes
nn_s = ne_x + 1;                                        % Number of nodes per side
dofn = 3;                                               % DOF per node
DOF = dofn*nn;                                          % DOF of the plate
structure.sets.Ug = (1:1:DOF)';
v_nn = (1:nn);
vdof = (1:DOF)';                                        % Vector containing all the DOF of the plate

% Load applied
P = -1000;                                          
xp = 0.5; yp = 0.5;                                     % Position of the applied load

f = 2000;                       % Maximum frequency
vf = (1:2000);                  % Vector of frequencies

% Location of the load
ixf = round(xp/plate.a/(plate.a/ne_x)) + 1;
iyf = round(yp/plate.b/(plate.b/ne_x)) + 1;
node_p = (ne_x + 1)*(iyf - 1) + ixf;

% Boundary conditions for the plate: ff (free-free) or ss (simply supported)
solve = 'ss';                

% Mesh the plate
[coordinates, nodes] = MeshRectanglularPlate(plate.a, plate.b,ne_x,ne_x);
PlotMesh(coordinates, nodes)

%% STIFFNESS, INERTIA AND FORCE MATRICES

coord = zeros(nn_s,2);                        % Coordinates of the external nodes
coord(:,1) = 0:plate.a/(nn_s-1):plate.a;      % X coordinates
coord(:,2) = 0:plate.b/(nn_s-1):plate.b;      % Y coordinates

[X,Y] = meshgrid(coord(:,1),coord(:,2));      % Mesh to represent the plate
Z = zeros(nn_s);
surf(X,Y,Z)

coord_n = zeros(nn_s,2*nn_s);
coord_n(:,1:2:end) = X;
coord_n(:,2:2:end) = Y;

coord_n = reshape(coord_n,[nn_s,2,nn_s]);       % Y coordinates are in the first column and X coordinates in the second

% Vector with the relative position of every element of the beam
row_col = zeros(ne,2);
for e  = 1:ne
    row_col(e,1) = floor((e-1)/ne_x) + 1;
    row_col(e,2) = mod((e-1),ne_x) + 1;
end

% Connectivity matrix
connectivity = zeros(2,2,ne);
for e = 1:ne
    node(1) = (ne_x + 1)*(row_col(e,1) - 1) + row_col(e,2);
    node(2) = (ne_x + 1)*(row_col(e,1) - 1) + row_col(e,2) + 1;
    node(3) = (ne_x + 1)*row_col(e,1) + row_col(e,2);
    node(4) = (ne_x + 1)*row_col(e,1) + row_col(e,2) + 1;

    connectivity(:,:,e) = [node(1) node(2); node(3) node(4)];
end

% Stiffness Matrix and Force vector
K = zeros(DOF); M = zeros(DOF);F = zeros(DOF,1);
D_b = plate.I*plate.E/(1 - plate.nu^2)*[1 plate.nu 0; plate.nu 1 0; 0 0 (1 - plate.nu)/2];  % The inertia is included to take into account the thikness of the plate
D_s = plate.G*plate.t*5/6*[1 0; 0 1];
index = zeros(2);

for e = 1:ne
   index = connectivity(:,:,e);
   x2 = coord_n(mod((e-1),ne_x) + 2,2,floor((e-1)/ne_x) + 1); x1 = coord_n(mod((e-1),ne_x) + 1,2,floor((e-1)/ne_x) + 1);
   y2 = coord_n(mod((e-1),ne_y) + 1,1,floor((e-1)/ne_y) + 2); y1 = coord_n(mod((e-1),ne_y) + 1,1,floor((e-1)/ne_y) + 1);
   xe = x2 - x1; ye = y2 - y1;                                                  % Lengths of the element on both directions
   detJe = xe*ye/4;                                                             % Transformation's Jacobian
   dofe = [index(1,1)*dofn - 2 index(1,1)*dofn - 1 index(1,1)*dofn...           % DOF of the element
           index(1,2)*dofn - 2 index(1,2)*dofn - 1 index(1,2)*dofn...
           index(2,2)*dofn - 2 index(2,2)*dofn - 1 index(2,2)*dofn...
           index(2,1)*dofn - 2 index(2,1)*dofn - 1 index(2,1)*dofn];

    % Bending matrix
    Kbe = zeros(dofn*nne);                                          % Initialization of element bending matrix

    for i = 1:iob
        xi = chi_ip.s(i);
        for j = 1:iob
            eta = chi_ip.s(j);                                      % Value of interpolation points
            weigth = w_ip.s(i);                                     % Weight of interpolation
            [~,dNdxi,dNdeta] = ShapeFunctions(xi,eta);              % Shape functions and derivatives 
            % Bending kinematic matrix
            B_b = [0 0 -dNdxi(1)*2/xe 0 0 -dNdxi(2)*2/xe 0 0 -dNdxi(3)*2/xe 0 0 -dNdxi(4)*2/xe;
                   0 dNdeta(1)*2/ye 0 0 dNdeta(2)*2/ye 0 0 dNdeta(3)*2/ye 0 0 dNdeta(4)*2/ye 0; 
                   0 dNdxi(1)*2/ye -dNdeta(1)*2/xe 0 dNdxi(2)*2/ye -dNdeta(2)*2/xe 0 dNdxi(3)*2/ye -dNdeta(3)*2/xe 0 dNdxi(4)*2/ye -dNdeta(4)*2/xe]; 
            % Bending matrix of the element
            Kbe = Kbe + B_b'*D_b*B_b*detJe*weigth*weigth;                        % The weights are in both xi and eta
        end
    end

    % Shear matrix
    Kse = zeros(dofn*nne);                              % Initialization of element shear matrix

    % Selective integration allows to integrate the shear component with order 1 instead of 2
    xi = chi_ip.f; eta = chi_ip.f;weigth = w_ip.f;
    [N,dNdxi,dNdeta] = ShapeFunctions(xi,eta);              % Shape functions and derivatives
    % Shear kinematic matrix
    B_s =[dNdxi(1)*2/xe 0 N(1) dNdxi(2)*2/xe 0 N(2) dNdxi(3)*2/xe 0 N(3) dNdxi(4)*2/xe 0 N(4);
          dNdeta(1)*2/ye -N(1) 0 dNdeta(2)*2/ye -N(2) 0 dNdeta(3)*2/ye -N(3) 0 dNdeta(4)*2/ye -N(4) 0];

    % Shear matrix of the element
    Kse = Kse + B_s'*D_s*B_s*detJe*weigth*weigth;

    % Stiffness matrix of the element
    Ke = Kbe + Kse;                                     
    
    % Global stiffness matrix
    K(dofe,dofe) = K(dofe,dofe) + Ke;

    % Inertia matrix (lumped mass)
    me = plate.rho*plate.t*xe*ye;                           % Mass of the element
    
    Me = me/4*[1 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0;...
                    0 0 0 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0;...
                    0 0 0 0 0 0 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0;...
                    0 0 0 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0];

    M(dofe,dofe) = M(dofe,dofe) + Me;

%     % Force vector
%     Fe = zeros(dofn*nne,1);                             % Initialization of element force vector
% 
%     xi = chi_ip.f; eta = chi_ip.f;weigth = w_ip.f;
%     N_1 = (1 - xi)*(1 - eta)/4; N_2 = (1 + xi)*(1 - eta)/4;         
%     N_3 = (1 + xi)*(1 + eta)/4; N_4 = (1 - xi)*(1 + eta)/4;
%     
%     % Force vector of the element
%     fe = [N_1*P 0 0 N_2*P 0 0 N_3*P 0 0 N_4*P 0 0]';
%     Fe = Fe + fe*weigth*weigth*detJe;
% 
%     % Global force vector
%     F(dofe,1) = F(dofe,1) + Fe;
end

% Load vector
F(3*node_p - 2,1) = P;

%% STATIC PROBLEM ([K]*u = [F])

u = zeros(dofn*nn,1);

% Boundary conditions for simply suported beam
k1 = find(coordinates == 0); k2 = find(coordinates == plate.a);
k3 = find(coordinates > 0 & coordinates < 1);

bc_X = k1(1:length(k1)/2);                  % Nodes along the X-axis
bc_Y = k1(length(k1)/2+1:end) - nn;         % Nodes along the Y-axis
bc_xL = k2(1:length(k2)/2);                 % Nodes along X = L
bc_yL = k2(length(k2)/2+1:end) - nn;        % Nodes along Y = L

rdof_ss = unique([bc_X.*dofn-2; bc_Y.*dofn-2; bc_xL.*dofn-2; bc_yL.*dofn-2],"rows");   % Restricted DOF
fdof_ss = vdof;                             % Free DOF

for i = 1:length(rdof_ss)
    ip = rdof_ss(i);
    fdof_ss(ip) = vdof(ip) - rdof_ss(i);
end
fdof_ss = find(fdof_ss > 0);

% Boundary conditions for free free beam
rdof_ff = 0;
fdof_ff = (1:DOF)';

switch solve
    case 'ff'
    u = K\F;                           % Displacements in free DOF for free-free plate
    case 'ss'
    K_FF = K(fdof_ss,fdof_ss);         % Stiffness matrix in free DOF
    K_FR = K(fdof_ss,rdof_ss);
    M_FF = M(fdof_ss,fdof_ss);         % Mass matrix in free DOF
    F_F = F(fdof_ss,1);                % Force vector in free DOF
    
    u_F = K_FF\F_F;                    % Displacements in free DOF
    u(fdof_ss,1) = u_F;
    
    F_R = K_FR'*u_F;                   % Reactions in supports 
end           

%% PLOTS OF STATIC PROBLEM

% Deformed Shape
[w,thetax,thetay] = mytable(nn,u,DOF);
x = coordinates(:,1); y = coordinates(:,2);
figure(1)
plot3(x,y,w,'.')
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')

% Contour
figure(2)
PlotFieldonMesh(coordinates,nodes,w)

figure(3)
PlotFieldonDefoMesh(coordinates,nodes,w,w)

%% DYNAMIC SYSTEM (q0[[K] - Ω^2[M]] = p0)

D_FF = zeros(length(fdof_ss),length(fdof_ss),f);                    % Dynamic stiffness matrix of free DOF
q0 = zeros(DOF,f); q0_F = zeros(length(fdof_ss),f);

for i = 1:f
    D_FF(:,:,i) = (K_FF - (2*pi*i)^2*M_FF);
    q0_F(:,i) = D_FF(:,:,i)\F_F;
    q0(fdof_ss,i) = q0_F(:,i);
end

Q0 = squeeze(q0(node_p*3-2,:));

%% PLOTS OF DYNAMIC SYSTEM

% Amplitude vs Frequency
figure(4)
semilogy(vf,abs(Q0))
title("Amplitude Bode Diagram")

% Angular offset vs Frequency
figure(5)
plot(f,unwrap(angle(Q0)))
