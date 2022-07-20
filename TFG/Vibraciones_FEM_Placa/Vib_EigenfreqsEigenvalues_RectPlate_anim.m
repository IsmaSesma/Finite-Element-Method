%% **********************************************************************************
%      NATURAL FREQUENCIES AND MODE SHAPES OF A RECTANGULAR PLATE (Lx x Ly)
%                           CC-BY-NC-ND, M. Chimeno,  2022-07
% Code based on
% KSSV (2022). Plate Bending (https://www.mathworks.com/matlabcentral/fileexchange/32029-plate-bending)
%              MATLAB Central File Exchange. Retrieved April 29, 2022.
% In this modification:
%  - Setting custom plate and mesh parameters added in this modification
%  - Setting distributed force or point force
%  - Addition of the inertia matrix (lumped mass model)
%_____________________________________________________________________________________
clear all; close all; clc;
addpath('auxfunctions\')
%_____________________________________________________________________________________
disp('__________________________________________________________________________________');
disp(' ');
disp('     NATURAL FREQUENCIES AND MODE SHAPES OF A RECTANGULAR PLATE (Lx x Ly)');
disp('__________________________________________________________________________________');
disp(' ');
disp('1 - Definition of the plate:');
Lx = input('    Length along X axis (m)?    [Default 0.5 m]: ');
if isempty(Lx); Lx = 0.5; end
Ly = input('    Length along Y axis (m)?    [Default 1.0 m]: ');
if isempty(Ly); Ly = 1; end
t  = input('    Thickness (m)?            [Default 0.005 m]: ');
if isempty(t); t = 0.005; end
boundcond = input('    Boundary conditions on the four edges: simply supported (s) or clamped (c)?    [Default s]: ','s');
if isempty(boundcond), boundcond = 's'; end
switch boundcond
    case 's'
        typeBC = 'ss-ss-ss-ss';
    case 'c'
        typeBC = 'c-c-c-c' ;
end

disp('2 - Material Definition:');
E = input('    Modulus of elasticity (GPa)?    [Default 70 GPa]: ');
if isempty(E); E = 70e9; else E = E*1e9; end
nu= input('    Poison ratio?                   [Default 0.3]: ');
if isempty(nu); nu = 0.3; end
rho= input('    Density?                        [Default 2750 kg/m^3]: ');
if isempty(rho); rho = 2750; end
% ____________________________________________________________________________________
%               Model System definition (Mesh size)
if Lx<=Ly
    Nx       = 10;       % number of elements along the longer side
    Ny       = round(Nx/Lx*Ly);
else
    Ny       = 10;       % number of elements along the longer side
    Nx       = round(Ny/Ly*Lx);
end
% -------------------------------------------------------------------------------------
I = t^3/12 ;
% ____________________________________________________________________________________
disp(' ');
disp('5 - Constructing FEM');
%% MESH DEFINITION
    structure = CreateMesh(Lx,Ly,Nx,Ny);
    nel = length(structure.mesh.elements.id) ;                  % number of elements
    nnel=4;                                % number of nodes per element
    ndof=3;                                % number of dofs per node
    nnode = length(structure.mesh.nodes.id) ;          % total number of nodes in system
    sdof=nnode*ndof;                       % total system dofs  
    edof=nnel*ndof;
    structure.sets.Ug = [1:1:sdof]';
%% STIFFNESS MATRIX 
    %--------------------------------------------------------------------------
    % Order of Gauss Quadrature
    %--------------------------------------------------------------------------
    nglb=1;                     % 2x2 Gauss-Legendre quadrature for bending 
    ngls=1;                     % 1x1 Gauss-Legendre quadrature for shear 
    %--------------------------------------------------------------------------
    % Initialization of matrices and vectors
    %--------------------------------------------------------------------------
    force = zeros(sdof,1) ;             % System Force Vector
    Kgg=zeros(sdof,sdof);         % system stiffness matrix
    index=zeros(edof,1);                % index vector
    B_pb=zeros(3,edof);                 % kinematic matrix for bending
    B_ps=zeros(2,edof);                 % kinematic matrix for shear
    %--------------------------------------------------------------------------
    %  Computation of element matrices and vectors and their assembly
    %--------------------------------------------------------------------------
    %  For bending stiffness
    [pointb,weightb] = GaussQuadrature('first');     % sampling points & weights
    D_pb= I*E/(1-nu*nu)*[1  nu 0; nu  1  0; 0  0  (1-nu)/2];  % bending material property
    %  For shear stiffness
    [points,weights] = GaussQuadrature('first');    % sampling points & weights
    G = 0.5*E/(1.0+nu);                             % shear modulus
    shcof = 5/6;
    shcof = 1;% shear correction factor
    D_ps=G*shcof*t*[1 0; 0 1];                      % shear material property
    % Calculating element matrices
    for iel=1:nel                        % loop for the total number of elements
        for i=1:nnel
            node(i)=structure.mesh.elements.nodes(iel,i);               % extract connected node for (iel)-th element
            xx(i)=structure.mesh.nodes.coords(node(i),1);       % extract x value of the node
            yy(i)=structure.mesh.nodes.coords(node(i),2);       % extract y value of the node
        end
        ke = zeros(edof,edof);              % initialization of element stiffness matrix 
        kb = zeros(edof,edof);              % initialization of bending matrix 
        ks = zeros(edof,edof);              % initialization of shear matrix 
        f = zeros(edof,1) ;                 % initialization of force vector                   
        %  Numerical integration for bending term
        for intx=1:nglb
            xi=pointb(intx,1);                     % sampling point in x-axis
            wtx=weightb(intx,1);                   % weight in x-axis
            for inty=1:nglb
                eta=pointb(inty,2);                    % sampling point in y-axis
                wty=weightb(inty,2) ;                  % weight in y-axis
                [shape,dhdr,dhds]=Shapefunctions(xi,eta);      % compute shape functions and derivatives at sampling point
                [detjacobian,invjacobian]=Jacobian(nnel,dhdr,dhds,xx,yy);  % compute Jacobian
                [dhdx,dhdy]=ShapefunctionDerivatives(nnel,dhdr,dhds,invjacobian);   % derivatives w.r.t. physical coordinate
                B_pb=PlateBending(nnel,dhdx,dhdy);    % bending kinematic matrix
                %  compute bending element matrix
                kb=kb+(B_pb'*D_pb*B_pb)*wtx*wty*detjacobian;
            end
        end 
        %  numerical integration for shear term
        for intx=1:ngls
            xi=points(intx,1);                  % sampling point in x-axis
            wtx=weights(intx,1);               % weight in x-axis
            for inty=1:ngls
                eta=points(inty,2);                  % sampling point in y-axis
                wty=weights(inty,2) ;              % weight in y-axis
                [shape,dhdr,dhds]=Shapefunctions(xi,eta);           % compute shape functions and derivatives at sampling point
                [detjacobian,invjacobian]=Jacobian(nnel,dhdr,dhds,xx,yy);  % compute Jacobian
                [dhdx,dhdy]=ShapefunctionDerivatives(nnel,dhdr,dhds,invjacobian);      % derivatives w.r.t. physical coordinate
                B_ps=PlateShear(nnel,dhdx,dhdy,shape);        % shear kinematic matrix
                %  compute shear element matrix
                ks=ks+(B_ps'*D_ps*B_ps)*wtx*wty*detjacobian;
            end
        end
        %  compute element matrix
        ke = kb+ks ;
        index=elementdof(node,nnel,ndof);% extract system dofs associated with element
        [Kgg]=assemblematrix(Kgg,ke,index);  % assemble element stiffness and force matrices 
    end

%% MASS MATRIX (LUMPED MASS MODEL)
    Mgg=zeros(sdof,sdof);         % system mass matrix
    for iel=1:nel                       % loop for the total number of elements
        for i=1:nnel
            node(i)=structure.mesh.elements.nodes(iel,i);               % extract connected node for (iel)-th element
            xx(i)=structure.mesh.nodes.coords(node(i),1);       % extract x value of the node
            yy(i)=structure.mesh.nodes.coords(node(i),2);       % extract y value of the node
        end
        Ae = [xx(2)-xx(1),yy(2)-yy(1)]*[xx(3)-xx(1),yy(3)-yy(1)]';
        m = rho*Ae*t;
        me = zeros(edof,edof);              % initialization of element stiffness matrix
        me(1,1) = m/4; me(4,4)=m/4; me(7,7)=m/4;me(10,10)=m/4; 
        index=elementdof(node,nnel,ndof);   % extract system dofs associated with element
        [Mgg]=assemblematrix(Mgg,me,index);  % assemble element stiffness and force matrices 
    end

    
%% APPLYING BOUNDARY CONDITIONS
    structure.sets.Us = BoundaryCondition(typeBC,structure.mesh.nodes.coords) ;
    structure.sets.Uf = setdiff(structure.sets.Ug,structure.sets.Us);
    Kff = Kgg(structure.sets.Uf,structure.sets.Uf);
    Mff = Mgg(structure.sets.Uf,structure.sets.Uf);

%% EIGENVALUES ANALYSIS
    [PSI,d] = eig(Kff,Mff);
    [f0,order] = sort(sqrt(sum(d,1))/2/pi());
    disp('Natural frequencies:')
    disp(['   ' num2str(f0(1)) ' Hz, ' num2str(f0(2)) ' Hz, ' num2str(f0(3)) ' Hz, ' num2str(f0(4)) ' Hz'])
    PSI_f = PSI(:,order);
    PSI_g = zeros(sdof,size(PSI_f,2));
    PSI_g(structure.sets.Uf,:) = PSI_f;

%% ************************************************************************
%                             PLOTTING
    %% EIGENMDOES (3D)
    figure('Color','white','units','normalized','outerposition',[0 0 1 1])
        modestoview = 4;
        Ts = 1./f0(1:modestoview);
        vt = [0:min(Ts)/20:1*min(Ts)];
        % Initialization of the required matrices
        for t=1:length(vt)
            for c=1:modestoview
                s(c) = subplot(2,modestoview/2,c);
                hold on
                out = PSI_g(1:3:end,c);
                maxout = max(max(abs(out)));
                X = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
                Y = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
                profile = zeros(size(structure.mesh.elements.nodes,2),length(structure.mesh.elements.nodes)) ;
                for iel=1:length(structure.mesh.elements.nodes)   
                     for i=1:size(structure.mesh.elements.nodes,2)
                        nd(i)=structure.mesh.elements.nodes(iel,i);       
                        X(i,iel)=structure.mesh.nodes.coords(nd(i),1);    
                        Y(i,iel)=structure.mesh.nodes.coords(nd(i),2);    
                     end   
                     profile(:,iel) = out(nd');         
                end
                profile = profile/maxout*0.1*max(Lx,Ly);
                fill3(X,Y,profile*sin(2*pi()*f0(c)*vt(t)),profile*sin(2*pi()*f0(c)*vt(t)))
%                 colorbar
                axis equal
                xlabel('x (m)')
                ylabel('y (m)')
                view(30,30)
                set(gca,'ZLim',[-0.1*max(Lx,Ly),0.1*max(Lx,Ly)])
                title(['Mode ' num2str(c) ', f_0=' num2str(round(f0(c),1)) ' Hz'])
                set(gca,'FontSize',14)
            end
            drawnow;
            if t<length(vt)
                for m = 1:modestoview
                    cla(s(m));
                end
            end
            
        end

     
   
        
        
    