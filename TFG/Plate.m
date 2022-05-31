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
plate.tmd = 0.01;                                    % Thickness of the TMD

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

%% INPUT DATA AND MESH

ne_x = 17;                                              % Number of elements in X
ne_y = round(ne_x/plate.a*plate.b);                    % Number of elements in Y (keep the elements as square as possible)
dofn = 3;                                               % DOF per node

cont = 1;
empty_elements = zeros(19,19); a = length(empty_elements)^2;
for i = 1:6     % Elements that are not in the model
    if i == 6
        empty = (i*ne_x+2):(i+1)*ne_x-1;
    else
        empty = (i*ne_x+2):2:((i+1)*ne_x-1);
    end
    empty_elements(i,1:length(empty)) = empty;
    
    j = i + 9;
    if j == ne_x-2
        empty = (j*ne_x+2):(j+1)*ne_x-1;
    else
        empty = (j*ne_x+2):2:((j+1)*ne_x-1);
    end
empty_elements(j,1:length(empty)) = empty;
end
    
empty_elements = reshape(empty_elements, [a,1]);
empty_elements = nonzeros(empty_elements);
empty_elements = 0;                                % Decomment to simulate uniform plate
tmd_elements = 0;                                  % Elements that act as Tunned Mass Dumper

structure = CreateMesh(plate.a, plate.b,ne_x,ne_y);                         % Mesh definition
PlotMesh(structure.mesh.nodes.coords,structure.mesh.elements.nodes)         % Mesh plot

nne = 4;                                                % Number of nodes per element
ne = length(structure.mesh.elements.id);                % Number of elements
nn = length(structure.mesh.nodes.id);                   % Number of nodes
nn_sx = ne_x + 1;                                       % Number of nodes per side in x
nn_sy = ne_y +1;                                        % Number of nodes per side in y 
DOF = dofn*nn;                                          % DOF of the plate
vdof = (1:DOF)';                                        % Vector containing all the DOF of the plate

% Load applied
P = -1;                                          
xp = 0; yp = 0;                                     % Position of the applied load (% of the length)

f = 2000;                                               % Maximum frequency
vf = (1:2000);                                          % Vector of frequencies

% Location of the load
ixf = round(xp*plate.a/(plate.a/ne_x)) + 1;
iyf = round(yp*plate.b/(plate.b/ne_y)) + 1;
node_p = (ne_x + 1)*(iyf - 1) + ixf;

% Boundary conditions for the plate: ff (free-free), cc (clamped) or ss (simply supported)
solve = 'ff';

%% CONNECTIVITY

coord_nx = 0:plate.a/(nn_sx-1):plate.a;         % X coordinates
coord_ny = 0:plate.b/(nn_sy-1):plate.b;         % Y coordinates
[X,Y] = meshgrid(coord_nx,coord_ny);            % Mesh of the plate

coord_n = zeros(nn_sy,2*nn_sx);
coord_n(:,1:2:end) = X;
coord_n(:,2:2:end) = Y;

coord_n = reshape(coord_n,[nn_sy,2,nn_sx]);       % X coordinates are in the first column and Y coordinates in the second

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

%% STIFFNESS, INERTIA AND FORCE MATRICES

% Initialization of required matrices
K = zeros(DOF); M = zeros(DOF);F = zeros(DOF,1);
D_b = plate.I*plate.E/(1 - plate.nu^2)*[1 plate.nu 0; plate.nu 1 0; 0 0 (1 - plate.nu)/2];  % The inertia is included to take into account the thikness of the plate
D_s = plate.G*plate.t*5/6*[1 0; 0 1];
index = zeros(2);

for e = 1:ne
   index = connectivity(:,:,e);

   x2 = coord_n(floor((e-1)/ne_x) + 2,1,mod((e-1),ne_x) + 2); x1 = coord_n(floor((e-1)/ne_x) + 2,1,mod((e-1),ne_x) + 1);
   y2 = coord_n(floor((e-1)/ne_x) + 2,2,mod((e-1),ne_x) + 2); y1 = coord_n(floor((e-1)/ne_x) + 1,2,mod((e-1),ne_x) + 2);
   
   xe = x2 - x1; ye = y2 - y1;                                                  % Lengths of the element on both directions
   detJe = xe*ye/4;                                                             % Transformation's Jacobian
   dofe = [index(1,1)*dofn - 2 index(1,1)*dofn - 1 index(1,1)*dofn...           % DOF of the element
           index(1,2)*dofn - 2 index(1,2)*dofn - 1 index(1,2)*dofn...
           index(2,2)*dofn - 2 index(2,2)*dofn - 1 index(2,2)*dofn...
           index(2,1)*dofn - 2 index(2,1)*dofn - 1 index(2,1)*dofn];

    if ismember(e,empty_elements)
        Kbe = zeros(dofn*nne);
        Kse = zeros(dofn*nne);
        Me = zeros(dofn*nne);
    else

        if ismember(e,tmd_elements)
            plate.t = plate.tmd;
            plate.rho = plate.m/plate.a/plate.b/plate.t;        % Plate's density
            plate.I = plate.t^3/12;
        else
            plate.t = plate.t;
        end

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

        % Inertia matrix (lumped mass)
        me = plate.rho*plate.t*xe*ye;                           % Mass of the element
        Me = zeros(dofn*nne);
        Me(1,1) = me/4; Me (4,4) = me/4; Me(7,7) = me/4; Me(10,10) = me/4;
    end

        % Stiffness matrix of the element
        Ke = Kbe + Kse;                                     
        % Global stiffness matrix
        K(dofe,dofe) = K(dofe,dofe) + Ke;
        % Global inertia matrix
        M(dofe,dofe) = M(dofe,dofe) + Me;
end

% Load vector
F(3*node_p - 2,1) = P;

%% STATIC PROBLEM ([K]*u = [F])

u = zeros(dofn*nn,1);

% Boundary conditions for simply suported beam
k1 = find(structure.mesh.nodes.coords == 0); k2 = find(structure.mesh.nodes.coords == plate.a);

bc_X = k1(1:nn_sy);                             % Nodes along X = 0
bc_Y = k1(nn_sy+1:end) - nn;                    % Nodes along Y = 0
bc_xL = k2(1:nn_sy);                            % Nodes along X = L
bc_yL = k2(nn_sy+1:end) - nn/nn_sy*nn_sx;       % Nodes along Y = L

switch solve
    case 'ff'
        % Boundary conditions
        rdof = 0;
        fdof = (1:DOF)';
        
    case 'ss'
        % Boundary conditions
        rdof = unique([bc_X.*dofn-2; bc_Y.*dofn-2; bc_xL.*dofn-2; bc_yL.*dofn-2],"rows");   % Restricted DOF
        fdof = vdof;                             % Free DOF
    
        for i = 1:length(rdof)
            ip = rdof(i);
            fdof(ip) = vdof(ip) - rdof(i);
        end
        fdof = find(fdof > 0);

    case 'cc'
        % Boundary conditios
        rdof = unique([bc_X.*dofn-2; bc_X.*dofn-1; bc_X.*dofn; ...
                      bc_Y.*dofn-2; bc_Y.*dofn-1; bc_Y.*dofn; ...
                      bc_xL.*dofn-2; bc_xL.*dofn-1; bc_xL.*dofn; ...
                      bc_yL.*dofn-2; bc_yL.*dofn-1; bc_yL.*dofn],"rows");
        fdof = vdof;                           
    
        for i = 1:length(rdof)
            ip = rdof(i);
            fdof(ip) = vdof(ip) - rdof(i);
        end
        fdof = find(fdof > 0);

end           

 % Solve the system
 K_FF = K(fdof,fdof);         % Stiffness matrix in free DOF
 %K_FR = K(fdof,rdof);
 M_FF = M(fdof,fdof);         % Mass matrix in free DOF
 F_F = F(fdof,1);                % Force vector in free DOF
     
 u_F = K_FF\F_F;                    % Displacements in free DOF
 u(fdof,1) = u_F;
        
 %F_R = K_FR'*u_F;                   % Reactions in supports 

w = u(1:3:DOF);
thetax = u(2:3:DOF);
thetay = u(3:3:DOF);

%% PLOTS OF STATIC PROBLEM

% Plots the displacement and angles of the plate in 2D
figure('Color','white','units','normalized','outerposition',[0.3 0.3 0.5 0.6])
t = tiledlayout(1,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';
labels= {'Vertical Displacement (mm)','\theta_x (º)','\theta_y (º)'};
% Initialization of the required matrices
out = [w,rad2deg(thetax),rad2deg(thetay)];
nd = zeros();
for c=1:3
    nexttile(c)
    X = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
    Y = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
    profile = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
    for iel=1:length(structure.mesh.elements.nodes)   
        for i=1:size(structure.mesh.elements.nodes,2)
        nd(i)=structure.mesh.elements.nodes(iel,i);       
        X(i,iel)=structure.mesh.nodes.coords(nd(i),1);    
        Y(i,iel)=structure.mesh.nodes.coords(nd(i),2);    
        end   
        profile(:,iel) = out(nd',c) ;         
    end
    fill(X,Y,profile)
    colorbar
    axis equal
    xlabel('x (m)')
    ylabel('y (m)')
    title(labels{c})
end

% Plots the 3D plate deformed
figure('Color','white','units','normalized','outerposition',[0.3 0.3 0.5 0.6])
X = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
Y = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
Z = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
profile = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
for iel=1:length(structure.mesh.elements.nodes)   
    for i=1:size(structure.mesh.elements.nodes,2)
    nd(i)=structure.mesh.elements.nodes(iel,i);       
    X(i,iel)=structure.mesh.nodes.coords(nd(i),1);    
    Y(i,iel)=structure.mesh.nodes.coords(nd(i),2);    
    end
    Z(:,iel)=w(nd);                 
end
fill3(X,Y,Z*1000,Z*1000)
colorbar
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (mm)')
title('Vertical displacement on the deformed mesh (mm)')

%%
figure()
plot3(structure.mesh.nodes.coords(:,1),structure.mesh.nodes.coords(:,2),w,'.')

%% DYNAMIC SYSTEM (q0[[K] - Ω^2[M]] = p0)

%D_FF = zeros(length(fdof),length(fdof),f);                    % Dynamic stiffness matrix of free DOF
q0 = zeros(DOF,f); q0_F = zeros(length(fdof),f);

for i = 1:f
    D_FF(:,:) = (K_FF - (2*pi*i)^2*M_FF);
    q0_F(:,i) = D_FF(:,:)\F_F;
    q0(fdof,i) = q0_F(:,i);
end

Q0 = squeeze(q0(node_p*3-2,:));

%% PLOTS OF DYNAMIC SYSTEM

% Amplitude vs Frequency
figure()
semilogy(vf,abs(Q0))
title("Amplitude Bode Diagram")

% Angular offset vs Frequency
figure()
plot(vf,unwrap(angle(Q0)))
