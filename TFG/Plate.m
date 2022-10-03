% ***********************************************************************
%           FEM MODEL OF A PLATE USING RECTANGULAR ELEMENTS 
%                 Ismael Rodríguez Sesma, ETSIAE
% ***********************************************************************

clc;clear;close all;

%% MASIC AND GEOMETRIC DATA

plate.E = 3E9; plate.nu = 0.3;
plate.a = 0.169; plate.b = 0.168;                             % Plate's dimensions
plate.t = 0.004;                                            % Plate's thickness
plate.m = 0.130;                                            % Plate's mass
plate.rho = plate.m/plate.a/plate.b/plate.t;                % Plate's density
plate.rhos = plate.rho*plate.t;                             % Plate's surface density
plate.rhol = plate.rhos*plate.b;
plate.I = plate.t^3/12;                                     % Plate's inertia
plate.G = plate.E/2/(1 + plate.nu);
plate.D = plate.E*plate.t^3/12/(1 - plate.nu^2);
plate.beam_width = 0.001;                                   % Width of the tongue                               
plate.Leissa = plate.a^2*sqrt(plate.rhos/plate.D);          % Parameter defined in Leissa
acc_mass = 0.009;                                           % Mass of the accelerometers

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
% tmd.L = 0.005;                                               % Element length
% tmd.Leff = 6*tmd.L + 4*tmd.L/2;                              % Effective length
% tmd.ml = plate.rho*plate.beam_width*(6*tmd.L)*(2*tmd.L);     % Mass of the tongue
% tmd.I = tmd.L*plate.beam_width^3/12;
% tmd.t = 0.004:0.001:0.1;                
% 
% w0 = zeros(length(tmd.t),1);
% for i = 1:length(tmd.t)
%      w0(i,:) = sqrt(3*plate.E*tmd.I/tmd.Leff^3/(33*tmd.ml/140 + plate.rho*(tmd.L)^2*tmd.t(i)));
% end
% 
% figure()
% title('Election of tip width')
% xlabel('Width of the extra mass')
% ylabel('Resonance frequency')
% plot(tmd.t,w0)
% 
% plate.tmd = 0.1;

%% INPUT DATA AND MESH

ne_x = 34;                                              % Number of elements in X
ne_y = round(ne_x/plate.a*plate.b);                     % Number of elements in Y (keep the elements as square as possible)
dofn = 3;                                               % DOF per node

% Mesh to mitigate first mode of the free-free beam
empty_elements = zeros(19,19);
for i = 1:11                  % Elements that are not in the model
    if i == 11
        empty = (i*ne_x+2):((i+1)*ne_x-5);
        no_empty(:,1) = (i*ne_x+6):6:((i+1)*ne_x-9);
        no_empty(:,2) = (i*ne_x+7):6:((i+1)*ne_x-9);

        empty = setxor(empty,no_empty);
    else
        empty = (i*ne_x+2):3:((i+1)*ne_x-5);
    end
    empty_elements(i,1:length(empty)) = empty;
    
    j = i + 21;               % Second row of TMD      
    if j == ne_x-2
        empty = (j*ne_x+2):((j+1)*ne_x-5);
        no_empty(:,1) = (j*ne_x+6):6:((j+1)*ne_x-9);
        no_empty(:,2) = (j*ne_x+7):6:((j+1)*ne_x-9);

        empty = setxor(empty,no_empty);
    else
        empty = (j*ne_x+2):3:((j+1)*ne_x-5);
    end
    empty_elements(j,1:length(empty)) = empty;
end

% Elements of the beams
beam_elements = zeros(); 
for i = 1:10               
    beam = (i*ne_x+3):((i+1)*ne_x-6);
    empty = (i*ne_x+5):3:((i+1)*ne_x-8);
    no_beam(:,1) = (i*ne_x+6):6:((i+1)*ne_x-10);
    no_beam(:,2) = (i*ne_x+7):6:((i+1)*ne_x-9);
 
    beam = setxor(beam,no_beam);
    beam = setxor(beam,empty);
    beam_elements(i,1:length(beam)) = beam;
    
    j = i + 21;                
    beam = (j*ne_x+3):((j+1)*ne_x-6);
    empty = (j*ne_x+5):3:((j+1)*ne_x-8);
    no_beam(:,1) = (j*ne_x+6):6:((j+1)*ne_x-10);
    no_beam(:,2) = ((j*ne_x+7):6:((j+1)*ne_x-9));

    beam = setxor(beam,no_beam);
    beam = setxor(beam,empty);
    beam_elements(j,1:length(beam)) = beam;
end

empty_elements = nonzeros(empty_elements);
beam_elements = nonzeros(beam_elements);

% Elements that act as Tunned Mass Dumper
% tmd_elements = [beam_elements(4:5:69),beam_elements(5:5:70)];
% tmd_elements = reshape(tmd_elements,[28,1]);

% Decomment to simulate uniform plate
empty_elements = 0; beam_elements = 0; tmd_elements = 0;                

% beam_elements = setxor(beam_elements,tmd_elements);
% acc_elements = [1 1136 1137 1138];                          % Elements where the accelerometers are placed
thread_nodes = [1124 1155];                                % Nodes where the thread is placed
thread_DOF = 3*thread_nodes;
acc_elements = [0 0 0 0];
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
P = -1;                                          
xp = 0; yp = 0;                                     % Position of the applied load (% of the length)

f = 2000;                                               % Maximum frequency
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

        if ismember(e,beam_elements)                                     % Recalculate Di for the elements with less mass
            plate.t = plate.beam;
            plate.I = plate.t^3/12;
            D_b = plate.I*plate.E/(1 - plate.nu^2)*[1 plate.nu 0; plate.nu 1 0; 0 0 (1 - plate.nu)/2];
            D_s = plate.G*plate.t*5/6*[1 0; 0 1];

        elseif ismember(e,tmd_elements)                                     % Recalculate Di for the elements with extra mass
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

        % Inertia matrix (lumped mass)
        me = plate.rho*plate.t*xe*ye;                           % Mass of the element
        Me = zeros(dofn*nne);
%         Me(1,1) = me/4; Me (4,4) = me/4; Me(7,7) = me/4; Me(10,10) = me/4;
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

%  % Solve the system
 K_FF = K(fdof,fdof);         % Stiffness matrix in free DOF
 K_FF_testK_FF(thread_DOF(1), thread_DOF(1)) = K_FF(thread_DOF(1), thread_DOF(1)) + 46E3;
 K_FF_testK_FF(thread_DOF(2), thread_DOF(2)) = K_FF(thread_DOF(2), thread_DOF(2)) + 46E3;
 %K_FR = K(fdof,rdof);
 M_FF = M(fdof,fdof);         % Mass matrix in free DOF
 F_F = F(fdof,1);                % Force vector in free DOF

 fprintf('Matrices computed\n');
     
%  u_F = K_FF\F_F;                    % Displacements in free DOF
%  u(fdof,1) = u_F;
%         
%  %F_R = K_FR'*u_F;                   % Reactions in supports 
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
% 
% q0 = zeros(DOF,f); q0_F = zeros(length(fdof),f);
% 
% for i = 1:f
%     D_FF(:,:) = (K_FF - (2*pi*i)^2*M_FF);
%     q0_F(:,i) = D_FF(:,:)\F_F;
%     q0(fdof,i) = q0_F(:,i);
% end

%% PLOTS OF DYNAMIC SYSTEM
% 
% % Amplitude vs Frequency
% Q0 = squeeze(q0(943,:));
% Q1 = squeeze(q0(1,:));
% figure()
% semilogy(vf,abs(Q0))
% title("Amplitude Bode Diagram")
% 
% figure()
% semilogy(vf,abs(Q1))
% 
% figure()
% semilogy(vf,abs(Q0./Q1))
% 
% % Angular offset vs Frequency
% figure()
% plot(vf,unwrap(angle(Q0)))

[PSI,d] = eig(K_FF,M_FF);
[f0,order] = sort(sqrt(sum(d,1))/2./pi);
disp('Natural frequencies:')
disp(['   ' num2str(f0(5)) ' Hz, ' num2str(f0(6)) ' Hz, ' num2str(f0(7)) ' Hz, ' num2str(f0(8)) ' Hz'])
PSI_f = PSI(:,order);
PSI_g = zeros(DOF,size(PSI_f,2));
PSI_g(fdof,:) = PSI_f;

fprintf('Natural frequencies and eigenmodes obtained\n');

%% PLOTS OF SHAPE MODES

% figure('Color','white','units','normalized','outerposition',[0 0 1 1])
% modestoview = 8;
% Ts = 1./f0(5:modestoview+4);
% vt = (0:min(Ts)/20:1*min(Ts)); 
% Initialization of the required matrices
% for t=1:length(vt)
%     for c=5:modestoview+4
%         s(c-4) = subplot(2,modestoview/2,c-4);
%         hold on
%         out = PSI_g(1:3:end,c);
%         maxout = max(max(abs(out)));
%         X = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
%         Y = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
%         profile = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
%             for iel=1:length(structure.mesh.elements.nodes)   
%                 for i=1:size(structure.mesh.elements.nodes,2)
%                     nd(i)=structure.mesh.elements.nodes(iel,i);       
%                     X(i,iel)=structure.mesh.nodes.coords(nd(i),1);    
%                     Y(i,iel)=structure.mesh.nodes.coords(nd(i),2);    
%                 end   
%                 profile(:,iel) = out(nd');         
%             end
%             profile = profile/maxout*0.1*max(plate.a,plate.b);
%             fill3(X,Y,profile*sin(2*pi()*f0(c)*vt(t)),profile*sin(2*pi()*f0(c)*vt(t)))
%           colorbar
%             axis equal
%             xlabel('x (m)')
%             ylabel('y (m)')
%             view(30,30)
%             set(gca,'ZLim',[-0.1*max(plate.a,plate.b),0.1*max(plate.a,plate.b)])
%             title(['Mode ' num2str(c) ', f_0=' num2str(round(f0(c),1)) ' Hz'])
%             set(gca,'FontSize',14)
%     end
%     drawnow;
%     if t<length(vt)
%         for m = 1:modestoview
%             cla(s(m));
%         end
%     end            
% end
% 
%% VALIDATION (LEISSA)

table_ss = [19.7329 49.348 49.348 78.9568];        % Leissa values for a ss plate
freq_ss = table_ss/plate.Leissa/2/pi;

table_ff = [13.489 19.789 24.432 35.024];          % Leissa values for a ff plate
freq_ff = table_ff/plate.Leissa/2/pi;

disp('Leissa frequencies for FFFF plate:')
disp(['    ' num2str(freq_ff(1)) ' Hz, '  num2str(freq_ff(2)) ' Hz  ', num2str(freq_ff(3)) ' Hz  ', num2str(freq_ff(4)) ' Hz  '])
