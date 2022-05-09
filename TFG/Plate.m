% ***********************************************************************
%           FEM MODEL OF A PLATE USING RECTANGULAR ELEMENTS 
% ***********************************************************************

% clc;clear;close all;

%% MASIC AND GEOMETRIC DATA

plate.E = 10920; plate.nu = 0.3;
plate.a = 1; plate.b = 1;                % Plate's dimensions
plate.t = 0.1;                              % Plate's thickness
plate.m = 0.5;                               % Plate's mass
plate.I = plate.t^3/12;                   % Plate's inertia
plate.G = plate.E/2/(1 + plate.nu);

%% INPUT DATA

ne_s = 20;                    % Number of bidimensional elements per side of the plate
nne = 4;                      % Number of nodes per element
nn_s = ne_s+1;                % Number of nodes per side
ne = ne_s^2;                  % Number of bidimensional elements in the whole plate
nn = nn_s^2;                  % Number of nodes
v_nn = (1:nn);
dofn = 3;                     % DOF per node
DOF = dofn*nn;                % DOF of the plate
vdof = (1:DOF)';              % Vector containing all the DOF of the plate

P = -1;                       % Uniform preassuer aplied on the plate
solve = 'ss';                 % Boundary conditions for the plate: ff (free-free) or ss (simply supported)

% Mesh the plate
[coordinates, nodes] = MeshRectanglularPlate(plate.a, plate.b,ne_s,ne_s);
PlotMesh(coordinates, nodes)

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
    row_col(e,1) = floor((e-1)/ne_s) + 1;
    row_col(e,2) = mod((e-1),ne_s) + 1;
end

% Connectivity matrix
connectivity = zeros(2,2,ne);
for e = 1:ne
    node1 = (ne_s + 1)*(row_col(e,1) - 1) + row_col(e,2);
    node2 = (ne_s + 1)*(row_col(e,1) - 1) + row_col(e,2) + 1;
    node3 = (ne_s + 1)*row_col(e,1) + row_col(e,2);
    node4 = (ne_s + 1)*row_col(e,1) + row_col(e,2) + 1;

    connectivity(:,:,e) = [node1 node2; node3 node4];
end

% Stiffness Matrix and Force vector
K = zeros(DOF); F = zeros(DOF,1); kk = zeros(DOF); ff = zeros(DOF,1);
D_b = plate.I*plate.E/(1 - plate.nu^2)*[1 plate.nu 0; plate.nu 1 0; 0 0 (1 - plate.nu)/2];  % The inertia is included to take into account the thikness of the plate
D_s = plate.G*plate.t*5/6*[1 0; 0 1];
index = zeros(2);

for e = 1:ne
   index = connectivity(:,:,e);
   x2 = coord_n(mod((e-1),ne_s) + 2,2,floor((e-1)/ne_s) + 1); x1 = coord_n(mod((e-1),ne_s) + 1,2,floor((e-1)/ne_s) + 1);
   y2 = coord_n(mod((e-1),ne_s) + 1,1,floor((e-1)/ne_s) + 2); y1 = coord_n(mod((e-1),ne_s) + 1,1,floor((e-1)/ne_s) + 1);
   xe = x2 - x1; ye = y2 - y1;                                                  % Lengths of the element on both directions
   detJe = xe*ye/4;                                                             % Transformation's Jacobian
   dofe = [index(1,1)*dofn - 2 index(1,1)*dofn - 1 index(1,1)*dofn...           % DOF of the element
           index(1,2)*dofn - 2 index(1,2)*dofn - 1 index(1,2)*dofn...
           index(2,1)*dofn - 2 index(2,1)*dofn - 1 index(2,1)*dofn...
           index(2,2)*dofn - 2 index(2,2)*dofn - 1 index(2,2)*dofn]; 

    % Bending matrix
    Kbe = zeros(dofn*nne);                                          % Initialization of element bending matrix

    for i = 1:iob
        xi = chi_ip.s(i);
        for ip = 1:iob
            eta = chi_ip.s(ip);                                      % Value of interpolation points
            weigth = w_ip.s(i);                                          % Weight of interpolation
            dNdxi_1 = -(1 - eta)/4; dNdeta_1 = -(1 - xi)/4;         % Derivates of shape functions
            dNdxi_2 = (1 - eta)/4; dNdeta_2 = -(1 + xi)/4; 
            dNdxi_3 = (1 + eta)/4; dNdeta_3 = (1 + xi)/4; 
            dNdxi_4 = -(1 + eta)/4; dNdeta_4 = (1 - xi)/4; 
            % Bending kinematic matrix
            B_b = [0 -dNdxi_1*2/xe 0 0 -dNdxi_2*2/xe 0 0 -dNdxi_3*2/xe 0 0 -dNdxi_4*2/xe 0;
                   0 0 -dNdeta_1*2/ye 0 0 -dNdeta_2*2/ye 0 0 -dNdeta_3*2/ye 0 0 -dNdeta_4*2/ye; 
                   0 -dNdeta_1*2/ye -dNdxi_1*2/xe 0 -dNdeta_2*2/ye -dNdxi_2*2/xe 0 -dNdeta_3*2/ye -dNdxi_3*2/xe 0 -dNdeta_4*2/ye -dNdxi_4*2/xe]; 
            % Bending matrix of the element
            Kbe = Kbe + B_b'*D_b*B_b*detJe*weigth*weigth;                        % The weights are in both xi and eta
        end
    end

    % Shear matrix
    Kse = zeros(dofn*nne);                              % Initialization of element shear matrix

    % Selective integration allows to integrate the shear component with order 1 instead of 2
    xi = chi_ip.f; eta = chi_ip.f;weigth = w_ip.f;

    N_1 = (1 - xi)*(1 - eta)/4; N_2 = (1 + xi)*(1 - eta)/4;           % Shape functions
    N_3 = (1 + xi)*(1 + eta)/4; N_4 = (1 - xi)*(1 + eta)/4;
    dNdxi_1 = -(1 - eta)/4; dNdeta_1 = -(1 - xi)/4;
    dNdxi_2 = (1 - eta)/4; dNdeta_2 = -(1 + xi)/4; 
    dNdxi_3 = (1 + eta)/4; dNdeta_3 = (1 + xi)/4; 
    dNdxi_4 = -(1 + eta)/4; dNdeta_4 = (1 - xi)/4;

    % Shear kinematic matrix
    B_s =[dNdxi_1*2/xe -N_1 0 dNdxi_2*2/xe -N_2 0 dNdxi_3*2/xe -N_3 0 dNdxi_4*2/xe -N_4 0;
          dNdeta_1*2/ye 0 -N_1 dNdeta_2*2/ye 0 -N_2 dNdeta_3*2/ye 0 -N_3 dNdeta_4*2/ye 0 -N_4];

    % Shear matrix of the element
    Kse = Kse + B_s'*D_s*B_s*detJe*weigth*weigth;

    Ke = Kbe + Kse;                                     % Stiffness matrix of the element
    
    % Global stiffness matrix
    K(dofe,dofe) = K(dofe,dofe) + Ke;

    % Force vector
    Fe = zeros(dofn*nne,1);                             % Initialization of element force vector

    xi = chi_ip.f; eta = chi_ip.f;weigth = w_ip.f;
    N_1 = (1 - xi)*(1 - eta)/4; N_2 = (1 + xi)*(1 - eta)/4;         
    N_3 = (1 + xi)*(1 + eta)/4; N_4 = (1 - xi)*(1 + eta)/4;
    
    % Force vector of the element
    fe = [N_1*P 0 0 N_2*P 0 0 N_3*P 0 0 N_4*P 0 0]';
    Fe = Fe + fe*weigth*weigth*detJe;

    % Global force vector
    F(dofe,1) = F(dofe,1) + Fe;
end

%% STATIC PROBLEM (K*u = F)

u = zeros(dofn*nn,1);

% Boundary conditions for simply suported beam
k1 = find(coordinates == 0); k2 = find(coordinates == plate.a);
k3 = find(coordinates > 0 & coordinates < 1);

bc_X = k1(1:length(k1)/2);                  % Nodes along the X-axis
bc_Y = k1(length(k1)/2+1:end) - nn;         % Nodes along the Y-axis
bc_xL = k2(1:length(k2)/2);                 % Nodes along X = L
bc_yL = k2(length(k2)/2+1:end) - nn;        % Nodes along Y = L

rdof_ss = unique([bc_X.*dofn-2; bc_X.*dofn; bc_Y.*dofn-2; bc_Y.*dofn-1; bc_xL.*dofn-2; bc_xL.*dofn; bc_yL.*dofn-2; bc_yL.*dofn-1],"rows");   % Restricted DOF
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
    F_F = F(fdof_ss,1);                % Force vector in free DOF
    
    u_F = K_FF\F_F;                    % Displacements in free DOF
    u(fdof_ss,1) = u_F;
    
    F_R = K_FR'*u_F;                   % Reactions in supports  
end           

%% PLOTS

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

%% FUNCTIONS

