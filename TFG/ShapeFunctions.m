function [shape,dNdxi,dNdeta] = ShapeFunctions(xi,eta)
%--------------------------------------------------------------------------
% Purpose:
%   Compute isoparametric shape functions of four node rectangular element
%   and their derivatives at a given integration point
%
% Variable Description:
%   shape: shape function for four-node element
%   dNdxi: derivatives of the shape function w.r.t. xi
%   dNdeta: derivatives of the shape function w.r.t. eta
%   xi: x coordinate value
%   eta: y coordinate value
%
%  Notes:
%     1st node at (-1,-1), 2nd node at (1,-1)
%     3rd node at (1,1), 4th node at (-1,1)
%--------------------------------------------------------------------------

% Shape functions
shape(1) = (1 - xi)*(1 - eta)/4;
shape(2) = (1 + xi)*(1 - eta)/4;
shape(3) = (1 + xi)*(1 + eta)/4;
shape(4) = (1 - xi)*(1 + eta)/4;

% Derivatives
dNdxi(1) = -(1 - eta)/4;
dNdxi(2) = (1 - eta)/4;
dNdxi(3) = (1 + eta)/4;
dNdxi(4) = -(1 + eta)/4;

dNdeta(1) = -(1 - xi)/4;
dNdeta(2) = -(1 + xi)/4;
dNdeta(3) = (1 + xi)/4;
dNdeta(4) = (1 - xi)/4;