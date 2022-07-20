function [kk]=assemblematrix(kk,k,index)
%----------------------------------------------------------
%  Purpose:
%     Assembly of element stiffness matrix into the system matrix &
%     Assembly of element force vector into the system vector
%
%  Synopsis:
%     [kk,ff]=assemble(kk,ff,k,index)
%
%  Variable Description:
%     kk - system stiffnes matrix
%     k  - element stiffness matrix
%     index - d.o.f. vector associated with an element
%-----------------------------------------------------------
edof = length(index);
for i=1:edof
    ii=index(i);
    for j=1:edof
        jj=index(j);
        kk(ii,jj)=kk(ii,jj)+k(i,j);
    end
end

