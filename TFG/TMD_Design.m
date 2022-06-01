% ***********************************************************
%           TMD DESIGN FOR THE PLATE ABSORBERS 
% ***********************************************************

% Provides with the dimensions of the dynamic absorber based on
% the frequency that wants to be eliminated

clc; clear; close all

%% INITIAL DATA

beam.E = 3E9;                                                   % Value based on experimental results for 3D printed beams
beam.dim.b = 0.01;                                              % Dimensions based on the plate model of 17 elements and 17x17 cm ===> The absorber is 1 element wide
beam.dim.t = 0.004; 
beam.I = beam.dim.b*beam.dim.t^3/12;

plate.dim.a = 0.17; plate.dim.b = 0.17; plate.dim.t = 0.004;
plate.m = 0.5; 
beam.rho = plate.m/plate.dim.a/plate.dim.b/plate.dim.t;         % Density if that of the plate   

beam.dim.L_abs = 0:0.001:0.05;                                  % Possible length of the absorber in the tip
beam.dim.L_beam = 0.05;                                         % Length of the beam part of the absorber
beam.dim.Leff = beam.dim.L_beam + beam.dim.L_abs/2;             % Possible effective lengths of the absorber
beam.dim.t_abs = 0.004:0.001:0.1;                             % Possible width of the absorber in the tip of the beam

%% NATURAL FREQUENCY (w0 = sqrt(Keq/Meq))

x = beam.dim.L_abs; y = beam.dim.t_abs;

w0 = zeros(length(y),length(x));
for i = 1:length(y)
    w0(i,:) = sqrt((3*beam.E*beam.I./(beam.dim.L_beam + x/2).^3)./(beam.rho*beam.dim.b*x*y(i)+33*beam.rho*beam.dim.t*beam.dim.L_beam/140));
end

figure(1)
surf(x,y,w0)
title('Natural frequency variation for L_{beam} = 5 cm')
xlabel('L_{extra mass}');ylabel('Width_{extra mass}');zlabel('\omega_0');

figure(2)
title('Natural frequency variation vs L_{extra mass}')
plot(x,w0)
xlabel('L_{extra mass}');ylabel('\omega_0');

figure(3)
title('Natural frequency variation vs Width_{extra mass}')
plot(y,w0)
xlabel('Width_{extra mass}');ylabel('\omega_0');

