function structure = CreateMesh(Lx,Ly,Nx,Ny)
    cont = 1;
    dof = 1;
    structure.mesh.nodes.id = [];
    cont = 1;
    for j=1:Ny+1
        for i=1:Nx+1
            structure.mesh.nodes.id(cont) = (Nx+1)*(j-1) + i;
            structure.mesh.nodes.coords(cont,:) = [(i-1)*(Lx/Nx),(j-1)*(Ly/Ny)];
            cont = cont + 1;
        end
    end
    cont = 1;
    for j=1:Ny
        for i=1:Nx
            structure.mesh.elements.id(cont) = (Nx+1)*(j-1) + i;
            n1 = (Nx+1)*(j-1) + i;
            n2 = (Nx+1)*(j-1) + i+1;
            n3 = (Nx+1)*(j+1-1) + i+1;
            n4 = (Nx+1)*(j+1-1) + i;
            structure.mesh.elements.nodes(cont,:) = [n1,n2,n3,n4] ;
            cont = cont + 1;
        end
    end