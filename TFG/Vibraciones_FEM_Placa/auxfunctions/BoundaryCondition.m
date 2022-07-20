function bcdof_out = BoundaryCondition(typeBC,coordinates)
% Code based on
% KSSV (2022). Plate Bending (https://www.mathworks.com/matlabcentral/fileexchange/32029-plate-bending)
%              MATLAB Central File Exchange. Retrieved April 29, 2022.
% In this modification:
%  - Modification to use structure variable and optimize code (code lines reduced to half)
%--------------------------------------------------------------------------
%   Purpose:
%           To determine the boundary conditions degree of freedom
%   Synopsis:
%           bcdof = BoundaryCondition(typeBC,coordinates)
%
%   Variable Description:
%           bcdof - boundary condition degree of freedom
%                   dof's are (UZ,RX,RY)
%           typeBC - string which gives the state of boundary condition
%           coordinates - geometric coordinates of nodes
%           
%--------------------------------------------------------------------------
bc.nodes{1} = find(coordinates(:,2)==min(coordinates(:,2))) ; % at y = 0 (along X-axes)
bc.nodes{2} = find(coordinates(:,2)==max(coordinates(:,2))) ; % at y = Ly (along X-axes)
bc.nodes{3} = find(coordinates(:,1)==min(coordinates(:,1))) ; % at x = 0 (along Y-axes)
bc.nodes{4} = find(coordinates(:,1)==max(coordinates(:,1))) ; % at x = Lx (along Y-axes)
switch typeBC
    case 'c-c-c-c'
        for l=1:4
            for n=1:length(bc.nodes{l})
                id = bc.nodes{l}(n);
                if n==1
                    bc.dof{l} = [3*id-2;3*id-1;3*id];
                else
                    bc.dof{l} = [bc.dof{l};3*id-2;3*id-1;3*id];
                end
            end
        end
    case 'ss-ss-ss-ss'
        for l=1:4
            for n=1:length(bc.nodes{l})
                id = bc.nodes{l}(n);
                if n==1
                    bc.dof{l} = [3*id-2];
                else
                    bc.dof{l} = [bc.dof{l};3*id-2];
                end
            end
        end

end
bcdof_out = union(bc.dof{1},bc.dof{2});
bcdof_out = union(bcdof_out,bc.dof{3});
bcdof_out = union(bcdof_out,bc.dof{4});