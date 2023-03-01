% ***********************************************************************
%           FEM MODEL OF A PLATE USING RECTANGULAR ELEMENTS 
%                 Ismael Rodríguez Sesma, ETSIAE
% ***********************************************************************

 clc;clear;close all;

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',2)

%% MASIC AND GEOMETRIC DATA

plate.E = 2.7E9; plate.nu = 0.3;
plate.a = 0.17; plate.b = 0.17;                           % Plate's dimensions
plate.t = 0.004;                                            % Plate's thickness
plate.rho = 1240;                % Plate's density
plate.m = plate.rho*plate.a*plate.b*plate.t;                                            % Plate's mass
plate.m = 0.137;
plate.rhos = plate.rho*plate.t;                             % Plate's surface density
plate.rhol = plate.rhos*plate.b;
plate.I = plate.t^3/12;                                     % Plate's inertia
plate.G = plate.E/2/(1 + plate.nu);
plate.D = plate.E*plate.t^3/12/(1 - plate.nu^2);                              
plate.Leissa = plate.a^2*sqrt(plate.rhos/plate.D);          % Parameter defined in Leissa
acc_mass = 0.009;                                           % Mass of the accelerometers
plate.thread_K = 2800;                                      % Stiffness of the thread in the test                                     

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

%% TMD DIMENSIONS AND DESIGN
% 
tmd.elements = 1; tmd.beam_elements = 4;                                            % Distribution of elements in the damper
plate.beam_thickness = 0.001;                                                       % Width of the tongue 
tmd.L = 0.005;                                                                      % Element length
tmd.width = 0.01;                                                                   % Element width
tmd.Leff = tmd.beam_elements*tmd.L + tmd.elements*tmd.L/2;                          % Effective length
tmd.ml = plate.rho*plate.beam_thickness*tmd.beam_elements*tmd.L*tmd.width;          % Mass of the tongue
tmd.I = tmd.width*plate.beam_thickness^3/12;                                        % Inertia of the tongue

tmd.t = 0.001:0.0001:0.1;                                                           % Thickness of the TMD (design variable)              
beam.t = 0.001:0.001:0.01;

w0 = zeros(length(tmd.t),length(beam.t),1);
for j = 1:length(beam.t)
    for i = 1:length(tmd.t)
        w0(i,j,:) = sqrt(3*plate.E*tmd.I/tmd.Leff^3/(33*plate.rho*beam.t(j)*tmd.beam_elements*tmd.L*tmd.width/140 + plate.rho*tmd.elements*tmd.L*tmd.width*tmd.t(i)))/2/pi;
    end
end

figure(1)
plot(tmd.t,w0(:,:,:))
set(gca,'XLim',[0, 0.02])
xlabel('Thicnkess of the extra mass [m]')

ylabel('Resonance frequency [Hz]')

plate.tmd = 0.0033;      % TMD thickness chosen

%% INPUT DATA AND MESH

ne_x = 34;                                              % Number of elements in X
ne_y = round(ne_x/plate.a*plate.b);                     % Number of elements in Y (keep the elements as square as possible)
dofn = 3;                                               % DOF per node

plate_type = 'homogeneous'; % Chose between locally resonant, holed or homogeneous for the meshing
switch plate_type
    case 'locally resonant'
        % Mesh to mitigate first mode of the free-free beam
        empty_elements = zeros(19,19);
        for i = 1:6                  % Elements that are not in the model
            if i == 6
                empty = ((i+1)*ne_x+4):((i+2)*ne_x-3);
                no_empty(:,1) = ((i+1)*ne_x+8):6:((i+2)*ne_x-8);
                no_empty(:,2) = ((i+1)*ne_x+9):6:((i+2)*ne_x-7);
        
                empty = setxor(empty,no_empty);
            else
                empty = ((i+1)*ne_x+4):3:((i+2)*ne_x-3);
            end
            empty_elements(i,1:length(empty)) = empty;
            
            j = i + 8;                % Second row of empty material      
            if j == 14
                empty = ((j+1)*ne_x+4):((j+2)*ne_x-3);
                no_empty(:,1) = ((j+1)*ne_x+8):6:((j+2)*ne_x-8);
                no_empty(:,2) = ((j+1)*ne_x+9):6:((j+2)*ne_x-7);
        
                empty = setxor(empty,no_empty);
            else
                empty = ((j+1)*ne_x+4):3:((j+2)*ne_x-3);
            end
            empty_elements(j,1:length(empty)) = empty;

            k = j + 8;                 % Third row of empty material 
            if k == 22
                empty = ((k+1)*ne_x+4):((k+2)*ne_x-3);
                no_empty(:,1) = ((k+1)*ne_x+8):6:((k+2)*ne_x-8);
                no_empty(:,2) = ((k+1)*ne_x+9):6:((k+2)*ne_x-7);
        
                empty = setxor(empty,no_empty);
            else
                empty = ((k+1)*ne_x+4):3:((k+2)*ne_x-3);
            end
            empty_elements(k,1:length(empty)) = empty;

            l = k + 8;                  % Fourth row of empty material
            if l == 30
                empty = ((l+1)*ne_x+4):((l+2)*ne_x-3);
                no_empty(:,1) = ((l+1)*ne_x+8):6:((l+2)*ne_x-8);
                no_empty(:,2) = ((l+1)*ne_x+9):6:((l+2)*ne_x-7);
        
                empty = setxor(empty,no_empty);
            else
                empty = ((l+1)*ne_x+4):3:((l+2)*ne_x-3);
            end

            empty_elements(l,1:length(empty)) = empty;
        end
        
        empty_elements = nonzeros(empty_elements);
        
        % Elements of the beams
        beam_elements = zeros(); 
        for i = 1:5                         % First row               
            beam = ((i+1)*ne_x+5):((i+2)*ne_x-4);
            empty = ((i+1)*ne_x+7):3:((i+2)*ne_x-6);
            no_beam(:,1) = ((i+1)*ne_x+8):6:((i+2)*ne_x-8);
            no_beam(:,2) = ((i+1)*ne_x+9):6:((i+2)*ne_x-7);
         
            beam = setxor(beam,no_beam);
            beam = setxor(beam,empty);
            beam_elements(i,1:length(beam)) = beam;
            
            j = i + 8;                     % Second row                
            beam = ((j+1)*ne_x+5):((j+2)*ne_x-4);
            empty = ((j+1)*ne_x+7):3:((j+2)*ne_x-6);
            no_beam(:,1) = ((j+1)*ne_x+8):6:((j+2)*ne_x-8);
            no_beam(:,2) = ((j+1)*ne_x+9):6:((j+2)*ne_x-7);
        
            beam = setxor(beam,no_beam);
            beam = setxor(beam,empty);
            beam_elements(j,1:length(beam)) = beam;

            k = j + 8;                     % Third row
            beam = ((k+1)*ne_x+5):((k+2)*ne_x-4);
            empty = ((k+1)*ne_x+7):3:((k+2)*ne_x-6);
            no_beam(:,1) = ((k+1)*ne_x+8):6:((k+2)*ne_x-8);
            no_beam(:,2) = ((k+1)*ne_x+9):6:((k+2)*ne_x-7);
        
            beam = setxor(beam,no_beam);
            beam = setxor(beam,empty);
            beam_elements(k,1:length(beam)) = beam;

            l = k + 8;                     % Fourth row
            beam = ((l+1)*ne_x+5):((l+2)*ne_x-4);
            empty = ((l+1)*ne_x+7):3:((l+2)*ne_x-6);
            no_beam(:,1) = ((l+1)*ne_x+8):6:((l+2)*ne_x-8);
            no_beam(:,2) = ((l+1)*ne_x+9):6:((l+2)*ne_x-7);
        
            beam = setxor(beam,no_beam);
            beam = setxor(beam,empty);
            beam_elements(l,1:length(beam)) = beam;
        end
        
        beam_elements = nonzeros(beam_elements);
        
        % Elements that act as Tunned Mass Dumper
        tmd_elements = [beam_elements(5:5:200)];

        beam_elements = setxor(beam_elements, tmd_elements);

    case 'holed'
        empty_elements = zeros(19,19);
        for i = 1:6                  % Elements that are not in the model
            empty = ((i+1)*ne_x+4):((i+2)*ne_x-3);
            no_empty(:,1) = ((i+1)*ne_x+8):6:((i+2)*ne_x-8);
            no_empty(:,2) = ((i+1)*ne_x+9):6:((i+2)*ne_x-7);
        
            empty = setxor(empty,no_empty);
            empty_elements(i,1:length(empty)) = empty;
            
            j = i + 8;                % Second row of empty material      
            
            empty = ((j+1)*ne_x+4):((j+2)*ne_x-3);
            no_empty(:,1) = ((j+1)*ne_x+8):6:((j+2)*ne_x-8);
            no_empty(:,2) = ((j+1)*ne_x+9):6:((j+2)*ne_x-7);
        
            empty = setxor(empty,no_empty);
            empty_elements(j,1:length(empty)) = empty;

            k = j + 8;                 % Third row of empty material 

            empty = ((k+1)*ne_x+4):((k+2)*ne_x-3);
            no_empty(:,1) = ((k+1)*ne_x+8):6:((k+2)*ne_x-8);
            no_empty(:,2) = ((k+1)*ne_x+9):6:((k+2)*ne_x-7);
        
            empty = setxor(empty,no_empty);
            empty_elements(k,1:length(empty)) = empty;

            l = k + 8;                  % Fourth row of empty material

            empty = ((l+1)*ne_x+4):((l+2)*ne_x-3);
            no_empty(:,1) = ((l+1)*ne_x+8):6:((l+2)*ne_x-8);
            no_empty(:,2) = ((l+1)*ne_x+9):6:((l+2)*ne_x-7);
        
            empty = setxor(empty,no_empty);
            empty_elements(l,1:length(empty)) = empty;
        end
        
        empty_elements = sort(nonzeros(empty_elements));

        dof_elements = zeros(19,19);          % Elements whose DOFs do exist (contact with plate)
        for i = 1:6                
            
            if i == 1
                empty = ((i+1)*ne_x+4):((i+2)*ne_x-3);
                no_empty(:,1) = ((i+1)*ne_x+8):6:((i+2)*ne_x-8);
                no_empty(:,2) = ((i+1)*ne_x+9):6:((i+2)*ne_x-7);
        
                empty = setxor(empty,no_empty);

            elseif i == 6
                empty = ((i+1)*ne_x+4):((i+2)*ne_x-3);
                no_empty(:,1) = ((i+1)*ne_x+8):6:((i+2)*ne_x-8);
                no_empty(:,2) = ((i+1)*ne_x+9):6:((i+2)*ne_x-7);
        
                empty = setxor(empty,no_empty);
            else
                empty = ((i+1)*ne_x+4):3:((i+2)*ne_x-3);
            end
            dof_elements(i,1:length(empty)) = empty;
            
            j = i + 8;                
            if j == 9
                empty = ((j+1)*ne_x+4):((j+2)*ne_x-3);
                no_empty(:,1) = ((j+1)*ne_x+8):6:((j+2)*ne_x-8);
                no_empty(:,2) = ((j+1)*ne_x+9):6:((j+2)*ne_x-7);
        
                empty = setxor(empty,no_empty);
            elseif j == 14
                empty = ((j+1)*ne_x+4):((j+2)*ne_x-3);
                no_empty(:,1) = ((j+1)*ne_x+8):6:((j+2)*ne_x-8);
                no_empty(:,2) = ((j+1)*ne_x+9):6:((j+2)*ne_x-7);
        
                empty = setxor(empty,no_empty);
            else
                empty = ((j+1)*ne_x+4):3:((j+2)*ne_x-3);
            end
            dof_elements(j,1:length(empty)) = empty;

            k = j + 8;
            if k == 17
                empty = ((k+1)*ne_x+4):((k+2)*ne_x-3);
                no_empty(:,1) = ((k+1)*ne_x+8):6:((k+2)*ne_x-8);
                no_empty(:,2) = ((k+1)*ne_x+9):6:((k+2)*ne_x-7);
        
                empty = setxor(empty,no_empty);
            elseif k == 22
                empty = ((k+1)*ne_x+4):((k+2)*ne_x-3);
                no_empty(:,1) = ((k+1)*ne_x+8):6:((k+2)*ne_x-8);
                no_empty(:,2) = ((k+1)*ne_x+9):6:((k+2)*ne_x-7);

                empty = setxor(empty,no_empty);
            else
                empty = ((k+1)*ne_x+4):3:((k+2)*ne_x-3);
            end
            dof_elements(k,1:length(empty)) = empty;

            l = k + 8;                 
            if l == 25
                empty = ((l+1)*ne_x+4):((l+2)*ne_x-3);
                no_empty(:,1) = ((l+1)*ne_x+8):6:((l+2)*ne_x-8);
                no_empty(:,2) = ((l+1)*ne_x+9):6:((l+2)*ne_x-7);
        
                empty = setxor(empty,no_empty);
            elseif l == 30
                empty = ((l+1)*ne_x+4):((l+2)*ne_x-3);
                no_empty(:,1) = ((l+1)*ne_x+8):6:((l+2)*ne_x-8);
                no_empty(:,2) = ((l+1)*ne_x+9):6:((l+2)*ne_x-7);
        
                empty = setxor(empty,no_empty);
            else
                empty = ((l+1)*ne_x+4):3:((l+2)*ne_x-3);
            end

            dof_elements(l,1:length(empty)) = empty;

        end
        
        dof_elements = nonzeros(dof_elements);
        nodof_elements = setxor(empty_elements, dof_elements);       % DOFs that do not exist (internal hole elements)
        beam_elements = 0; tmd_elements = 0;

    case 'homogeneous'
        empty_elements = 0; beam_elements = 0; 
        tmd_elements = 0;        
        nodof_elements = 0;
end

test = 'no';     % Chose between yes or no to simulate the aditional elements in the test
switch test

    case 'yes'
        acc_elements = [1 1136 1137 1138];                          % Elements where the accelerometers are placed
        thread_nodes = [1124 1155];                                 % Nodes where the thread is placed
        thread_DOF = 3*thread_nodes;

    case 'no'
        acc_elements = [0 0 0 0]; 
        thread_nodes = [0 0]; 
        thread_DOF = 3*thread_nodes;
end

% General mesh
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
P = 1;                                          
xp = 0; yp = 0;                                     % Position of the applied load (% of the length)

f = 1000;                                               % Maximum frequency
vf = (1:f);                                             % Vector of frequencies

% Location of the load
ixf = round(xp*plate.a/(plate.a/ne_x)) + 1;
iyf = round(yp*plate.b/(plate.b/ne_y)) + 1;
node_p = (ne_x + 1)*(iyf - 1) + ixf;

% Boundary conditions for the plate: ff (free-free), cc (clamped) or ss (simply supported)
solve = 'ff';

fprintf('Mesh finished\n');
 
%% CONNECTIVITY

coord_nx = 0:plate.a/(nn_sx-1):plate.a;         % X coordinates
coord_ny = 0:plate.b/(nn_sy-1):plate.b;         % Y coordinates
[X,Y] = meshgrid(coord_nx,coord_ny);            % Mesh of the plate

coord_n = zeros(nn_sy,2*nn_sx);
coord_n(:,1:2:end) = X;
coord_n(:,2:2:end) = Y;

coord_n = reshape(coord_n,[nn_sy,2,nn_sx]);       % X coordinates are in the first column and Y coordinates in the second

% Vector with the relative position of every element of the plate
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
D_b = plate.I*plate.E/(1 - plate.nu^2)*[1 plate.nu 0; plate.nu 1 0; 0 0 (1 - plate.nu)/2];  % The inertia is included to take into account the thicness of the plate
D_s = plate.G*plate.t*5/6*[1 0; 0 1];
index = zeros(2);

for e = 1:ne
   index = connectivity(:,:,e);

   x2 = coord_n(floor((e-1)/ne_x) + 2,1,mod((e-1),ne_x) + 2); x1 = coord_n(floor((e-1)/ne_x) + 2,1,mod((e-1),ne_x) + 1);
   y2 = coord_n(floor((e-1)/ne_x) + 2,2,mod((e-1),ne_x) + 2); y1 = coord_n(floor((e-1)/ne_x) + 1,2,mod((e-1),ne_x) + 2);
   
   xe = x2 - x1; ye = y2 - y1;                                                  % Length of the element on both directions
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

        if ismember(e,beam_elements)                                     % Recalculate Di for the elements with less thicnkess
            plate.t = plate.beam_thickness;
            plate.I = plate.t^3/12;
            D_b = plate.I*plate.E/(1 - plate.nu^2)*[1 plate.nu 0; plate.nu 1 0; 0 0 (1 - plate.nu)/2];
            D_s = plate.G*plate.t*5/6*[1 0; 0 1];

        elseif ismember(e,tmd_elements)                                  % Recalculate Di for the elements with extra thickness
            plate.t = plate.tmd;
            plate.I = plate.t^3/12;
            D_b = plate.I*plate.E/(1 - plate.nu^2)*[1 plate.nu 0; plate.nu 1 0; 0 0 (1 - plate.nu)/2];
            D_s = plate.G*plate.t*5/6*[1 0; 0 1];

        else
            plate.t = plate.t;
            plate.I = plate.t^3/12;
            D_b = plate.I*plate.E/(1 - plate.nu^2)*[1 plate.nu 0; plate.nu 1 0; 0 0 (1 - plate.nu)/2];
            D_s = plate.G*plate.t*5/6*[1 0; 0 1];
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

        % Inertia matrix (lumped mass matrix)
        me = plate.rho*plate.t*xe*ye;                           % Mass of the element
        Me = zeros(dofn*nne);

        if isequal(e,acc_elements(1))
            Me(1,1) = me/4 + acc_mass;

        elseif isequal(e,acc_elements(2))
            Me(10,10) = me/4 + acc_mass;
            
        elseif isequal(e,acc_elements(3))
            Me(7,7) = me/4 + acc_mass/2;
            Me(10,10) = me/4 + acc_mass/2;
            
        elseif isequal(e,acc_elements(4))
            Me(7,7) = me/4 + acc_mass;

        else
            Me(1,1) = me/4; Me (4,4) = me/4; Me(7,7) = me/4; Me(10,10) = me/4;

        end
    end

    % Stiffness matrix of the element
    Ke = Kbe + Kse;                     
    
    % Global stiffness matrix
    K(dofe,dofe) = K(dofe,dofe) + Ke;
    % Global inertia matrix
    M(dofe,dofe) = M(dofe,dofe) + Me;
end

% Eliminate DOFs of empty nodes (empty rows and columns)

switch plate_type

    case 'holed' 
        no_dof = zeros(length(nodof_elements), nne*dofn);

        for i = 1:length(nodof_elements)
            e_1 = nodof_elements(i);
            index = connectivity(:,:,e_1);

            dofe = [index(1,1)*dofn - 2 index(1,1)*dofn - 1 index(1,1)*dofn...
                index(1,2)*dofn - 2 index(1,2)*dofn - 1 index(1,2)*dofn...
                index(2,2)*dofn - 2 index(2,2)*dofn - 1 index(2,2)*dofn...
                index(2,1)*dofn - 2 index(2,1)*dofn - 1 index(2,1)*dofn];

            no_dof(i,:) = dofe;
        end

        no_dof = unique((reshape(no_dof,[],1)));

        % Reshape matrices

        K(no_dof,:) = [];
        K(:,no_dof) = [];

        M(no_dof,:) = [];
        M(:,no_dof) = [];

        DOF = size(K,1);           % New number of DOF

end

% Load vector
F(3*node_p - 2,1) = P;

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
 K_FF = K(fdof,fdof);         % Stiffness matrix for free DOF

 % Add the threads to simulate the real test
%  K_FF(thread_DOF(1) - 2, thread_DOF(1) - 2) = K_FF(thread_DOF(1) - 2, thread_DOF(1) - 2) + plate.thread_K/2;
%  K_FF(thread_DOF(2) - 2, thread_DOF(2) - 2) = K_FF(thread_DOF(2) - 2, thread_DOF(2) - 2) + plate.thread_K/2;

%  K_FR = K(fdof,rdof);         % If restricted DOF are needed
 M_FF = M(fdof,fdof);         % Mass matrix for free DOF
 F_F = F(fdof,1);             % Force vector for free DOF

 fprintf('Matrices computed\n');
     
%% Solve the system for restricted  DOF 
% u_F = K_FF\F_F;                    % Displacements in free DOF
% u(fdof,1) = u_F;
%         
% F_R = K_FR'*u_F;                   % Reactions in supports 
% 
% w = u(1:3:DOF);
% thetax = u(2:3:DOF);
% thetay = u(3:3:DOF);

%% PLOTS OF STATIC PROBLEM
% 
% % Plots the displacement and angles of the plate in 2D
% figure('Color','white','units','normalized','outerposition',[0.3 0.3 0.5 0.6])
% t = tiledlayout(1,3);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
% labels= {'Vertical Displacement (mm)','\theta_x (º)','\theta_y (º)'};
% % Initialization of the required matrices
% out = [w,rad2deg(thetax),rad2deg(thetay)];
% nd = zeros();
% for c=1:3
%     nexttile(c)
%     X = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
%     Y = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
%     profile = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
%     for iel=1:length(structure.mesh.elements.nodes)   
%         for i=1:size(structure.mesh.elements.nodes,2)
%         nd(i)=structure.mesh.elements.nodes(iel,i);       
%         X(i,iel)=structure.mesh.nodes.coords(nd(i),1);    
%         Y(i,iel)=structure.mesh.nodes.coords(nd(i),2);    
%         end   
%         profile(:,iel) = out(nd',c) ;         
%     end
%     fill(X,Y,profile)
%     colorbar
%     axis equal
%     xlabel('x (m)')
%     ylabel('y (m)')
%     title(labels{c})
% end
% 
% % Plots the 3D plate deformed
% figure('Color','white','units','normalized','outerposition',[0.3 0.3 0.5 0.6])
% X = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
% Y = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
% Z = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
% profile = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
% for iel=1:length(structure.mesh.elements.nodes)   
%     for i=1:size(structure.mesh.elements.nodes,2)
%     nd(i)=structure.mesh.elements.nodes(iel,i);       
%     X(i,iel)=structure.mesh.nodes.coords(nd(i),1);    
%     Y(i,iel)=structure.mesh.nodes.coords(nd(i),2);    
%     end
%     Z(:,iel)=w(nd);                 
% end
% fill3(X,Y,Z*1000,Z*1000)
% colorbar
% xlabel('x (m)')
% ylabel('y (m)')
% zlabel('z (mm)')
% title('Vertical displacement on the deformed mesh (mm)')
% 
% figure()
% plot3(structure.mesh.nodes.coords(:,1),structure.mesh.nodes.coords(:,2),w,'.')

%% DYNAMIC SYSTEM (q0[[K] - Ω^2[M]] = p0)

q0 = zeros(DOF,f); q0_F = zeros(length(fdof),f);
acc = zeros(length(fdof),f);
D_FF = zeros(DOF, DOF);

for i = 1:f
    D_FF(:,:) = (K_FF - (2*pi*i)^2*M_FF);
    q0_F(:,i) = D_FF(:,:)\F_F;
    q0(fdof,i) = q0_F(:,i);
    acc(fdof,i) = q0_F(:,i)/P*i^2;
end
fprintf('Dynamic system computed\n')
%%
%Q4 = squeeze(acc(1,:));
Q4 = squeeze(acc(1,:));

figure()
semilogy(vf,abs(Q4))
hold on
% semilogy(vf,abs(Q1))
% hold on
% semilogy(vf,abs(Q3))
hold on
semilogy(vf,abs(Q2))
set(gca,'XLim',[0, 400])
xlabel('Frequency [Hz]');
ylabel('Accelerance [$m/s^2/N$]')
legend( 'Homogeneous plate', 'Locally resonant plate for 250 Hz')

return 
%% PLOTS OF DYNAMIC SYSTEM
 
% Amplitude vs Frequency
Q0 = squeeze(q0(943,:));

% figure()
% semilogy(vf,abs(Q1))
% set(gca,'XLim',[0, 400])
% xlabel('Frequency [Hz]');
% ylabel('Accelerance [$m/s^2/N$]')

figure()
semilogy(vf,abs(Q2))
set(gca,'XLim',[0, 400])
xlabel('Frequency [Hz]');
ylabel('Accelerance [$m/s^2/N$]')

% figure()
% semilogy(vf,abs(Q0./Q1))

% Angular offset vs Frequency
% figure()
% plot(vf,unwrap(angle(Q0)))

fprintf('Plots of a displacement obtained\n')

% Solve the system via eigs

[PSI,d] = eigs(K_FF,M_FF, 8, "smallestabs");          % 8 eigs are computed to obtain first 4 non-zero
[f0,order] = sort(sqrt(sum(d,1))/2./pi);              % Sort in order and transform in Hz

disp('Natural frequencies:')
disp([num2str(f0(5)) ' Hz, ' num2str(f0(6)) ' Hz, ' num2str(f0(7)) ' Hz, ' num2str(f0(8)) ' Hz'])

% Eigenvectors
PSI_f = PSI(:,order);
PSI_g = zeros(DOF,size(PSI_f,2));
PSI_g(fdof,:) = PSI_f;                  % Shape modes                              

fprintf('Natural frequencies and eigenvectors obtained\n');

%% PLOTS OF SHAPE MODES

% 3D plot
% 
figure('Color','white','units','normalized','outerposition',[0 0 1 1])
modestoview = 6;
Ts = 1./f0(3:size(f0,2));
vt = (0:min(Ts)/20:1*min(Ts)); 

% for t=1:length(vt)                                  % Loop to animate the movement
    for c=5:modestoview+2                           % For each mode considered  
        s(c-2) = subplot(2,2,c-4);      % Plot dimension
        hold on
        out = PSI_g(1:3:end,c);                     % Components of mode shapes on vertical DOF
        maxout = max(max(abs(out)));
        
        X = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes));
        Y = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes));
        profile = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes));       % Each node, each element

        % Loop for each mode shape that defines the position of each 
        % point and for each point its value in the mode shape vector

        for iel=1:length(structure.mesh.elements.nodes)     % Each element is run 
            for i=1:size(structure.mesh.elements.nodes,2)       % Each node of the element is run
                nd(i)=structure.mesh.elements.nodes(iel,i);          % Node in consideration       
                X(i,iel)=structure.mesh.nodes.coords(nd(i),1);       % X-coordinate of the node       
                Y(i,iel)=structure.mesh.nodes.coords(nd(i),2);       % Y-coordinate of the node  
            end   
            profile(:,iel) = out(nd');      % Vector containing the value of the mode shape for each node         
        end

        profile = profile/maxout*0.1*max(plate.a,plate.b);      % Rescalate the vector
       %Decoment to see animation 
       %fill3(X,Y,profile*sin(2*pi*f0(c)*vt(t)),profile*sin(2*pi*f0(c)*vt(t)))
      
        fill3(X,Y,profile,profile*sin(2*pi*f0(c)))
        colorbar
        axis equal
        xlabel('x (m)',Interpreter='latex')
        ylabel('y (m)',Interpreter='latex')
        view(30,30)
        set(gca,'ZLim',[-0.1*max(plate.a,plate.b),0.1*max(plate.a,plate.b)])
        title(['$f_0=$' num2str(round(f0(c))) ' Hz'])
        set(gca,'FontSize',14,'TickLabelInterpreter','latex')
     end
%     drawnow;                            % Limited to 20 fps
%     if t<length(vt)
%         for m = 1:modestoview
%             cla(s(m));
%         end
%     end            
% end


%% 2D plot

figure('Color','white','units','normalized','outerposition',[0 0 1 1])
modestoview = 6;

for c=5:modestoview+2                           
    s(c-2) = subplot(2,modestoview/2,c-2);      
    hold on
    out = PSI_g(1:3:end,c);                     
    maxout = max(max(abs(out)));
        
    X = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes));
    Y = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes));
    profile = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)); 
    
    x_plot = linspace(0,plate.a,35);
    y_plot = linspace(0, plate.b, 35);

    [Xp,Yp] = meshgrid(x_plot, y_plot);
    Zp = zeros(35,35);
    
    for iel=1:length(structure.mesh.elements.nodes)     
        for i=1:size(structure.mesh.elements.nodes,2)   
            nd(i)=structure.mesh.elements.nodes(iel,i); 
            X(i,iel)=structure.mesh.nodes.coords(nd(i),1);
            Y(i,iel)=structure.mesh.nodes.coords(nd(i),2);
        end   
        profile(:,iel) = out(nd');     
    end

    %profile = profile/maxout*0.1*max(plate.a,plate.b);      % Rescalate the vector
    contour(X,Y,profile)
%   contour(X,Y,profile,[0 0 0],ShowText="on")
    axis equal
    xlabel('x (m)',Interpreter='latex')
    ylabel('y (m)',Interpreter='latex')
    title(['Mode ' num2str(c) ', $f_0=$' num2str(round(f0(c),1)) ' Hz'],Interpreter="latex")
    set(gca,'FontSize',14,'TickLabelInterpreter','latex')
end


%% VALIDATION (LEISSA)

table_ss = [19.7329 49.348 49.348 78.9568];        % Leissa values for a ss plate
freq_ss = table_ss/plate.Leissa/2/pi;

table_ff = [13.489 19.789 24.432 35.024];          % Leissa values for a ff plate
freq_ff = table_ff/plate.Leissa/2/pi;

disp('Leissa frequencies for FFFF plate:')
disp(['    ' num2str(freq_ff(1)) ' Hz, '  num2str(freq_ff(2)) ' Hz  ', num2str(freq_ff(3)) ' Hz  ', num2str(freq_ff(4)) ' Hz  '])

%% EXPERIMENTAL RESULTS

f = [134 254 324];                  % Resonance frequencies obtained in test
E_ref_a = E_ISO(plate.a,plate.rhol,plate.I,f);
E_ref_b = E_ISO(plate.b,plate.rhol,plate.I,f);

