% ***********************************************************************
%                   DATA PROCESSING OF DYNAMIC TESTING
% ***********************************************************************

% Collects all the numeric data from the beam's testing and places them
% into different structures to automate the post-processing

clc;clear;close all

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',2)

%% BEAM'S PROPERTIES

beam.input = importdata('geomdata.txt',';',2);               % CHANGE THE FILE
beam.m = beam.input.data(:,1);
beam.L = beam.input.data(:,2);
beam.b = beam.input.data(:,3);
beam.t = beam.input.data(:,4);
beam.Ixx = beam.b.*beam.t.^3./12;                      
beam.rho = beam.m./beam.L./beam.b./beam.t;             
beam.rhom = beam.rho.*beam.b.*beam.t;                 

%% CREATE TREE STRUCTURE

% 3-D printing parameters

beam.seed{1} = {'F','V','E'};                                       % Printing orientation (Flat, Vertical, On Edge)
beam.seed{2} = {'99','90' '85','70','60','50','40'};                % Infill percentage (%) 
beam.seed{3} = {'00','15','30','45','60','90'};                     % Raster angle (deg)
beam.seed{4} = {'1','2','3','4'};                                   % Layer thickness (1E-4 m)
beam.seed{5} = {'ZZ','T'};                                          % Infill pattern (Zig-Zag, Triangular)   
beam.seed{6} = {'a','b','c'};                                       % Sets of the same beam (3 models of each one)

 % Specify the files that are going to be read here
% casestoread = {'F40153Tb','F40153Tc','F40154Ta','F40154Tb','F40154Tc','F50903Ta','F50903Tb','F50903Tc','F50904Ta','F50904Tb','F50904Tc','F70603Ta','F70603Tb','F70603Tc','F85453ZZa','F85453ZZb','F85453ZZc','F99004ZZa','F99004ZZb','F99004ZZc'};
% casestoread = {'F40001Ta','F40001Tb','F40001Tc','F40002Ta','F40002Tb','F40002Tc','F50451Ta','F50451Tb','F50451Tc','F50902Ta','F50902Tb','F50902Tc','F70602Ta','F70602Tb','F70602Tc','F70604Ta','F70604Tb','F70604Tc','F85452ZZa','F85452ZZb','F85452ZZc','F85454ZZa','F85454ZZb','F85454ZZc','F90001ZZa','F90001ZZb','F90001ZZc','F90302ZZa','F90302ZZb','F90302ZZc','F90303ZZa','F90303ZZb','F90303ZZc','F90304ZZa','F90304ZZb','F90304ZZc','F99001ZZa','F99001ZZb','F99001ZZc','F99002ZZa','F99002ZZb','F99002ZZc','F99003ZZa','F99003ZZb','F99003ZZc'};
% casestoread = {'F90304ZZa','F90304ZZb','F90304ZZc'};
%casestoread = {'F40154Ta','F40154Ta','F40154Ta'};
casestoread = {'Ensayo_Pletina_Blanca'};

%% READ FILES

for s=1:length(casestoread)
    filecontent = readmatrix([casestoread{s} '.csv'],'Delimiter',';','NumHeaderLines',20,'DecimalSeparator',',');             % CHANGE THE FILE
    test.(casestoread{s}).freq = filecontent(:,1);
    test.(casestoread{s}).real = filecontent(:,2);
    test.(casestoread{s}).imag = filecontent(:,3);
    
    test.(casestoread{s}).amp = sqrt(filecontent(:,2).^2 + filecontent(:,3).^2);
    
    % Findpeaks

    tolook_f = test.(casestoread{s}).freq(10:end);
    tolook_v = log(test.(casestoread{s}).amp(10:end));
    [test.(casestoread{s}).peak, test.(casestoread{s}).rf] = findpeaks(tolook_v,tolook_f,'MinPeakDistance',800,'MinPeakProminence',.1,'NPeaks',3);
    tolook_f = test.(casestoread{s}).freq(10:end);
    tolook_v = log(1./test.(casestoread{s}).amp(10:end));
    [test.(casestoread{s}).min, test.(casestoread{s}).af] = findpeaks(tolook_v,tolook_f,'MinPeakDistance',500,'MinPeakProminence',.6,'NPeaks',3);                                           % Find the antiresonance frequencies
    test.(casestoread{s}).min = abs(test.(casestoread{s}).min);

    % Compute the elastic module with three different methods

    test.(casestoread{s}).E_iso = diag(E_ISO(beam.L(s),beam.rhom(s),beam.Ixx(s),test.(casestoread{s}).af));
    test.(casestoread{s}).E_af = diag(E_AFF(beam.L(s),beam.rhom(s),beam.Ixx(s),test.(casestoread{s}).af));
    test.(casestoread{s}).E_rf = diag(E_RFF(beam.L(s),beam.rhom(s),beam.Ixx(s),test.(casestoread{s}).rf));

end

%% FIGURES
figure()
tiledlayout(1,1)
for t=1:round(length(casestoread))
    nexttile
    hold on
    grid on
    plot(test.(casestoread{t}).freq,log(test.(casestoread{t}).amp))
%     plot(test.(casestoread{t}).rf,test.(casestoread{t}).peak,'ro')
%     plot(test.(casestoread{t}).af,test.(casestoread{t}).min,'bo')
    set(gca, 'XLim', [10,3000])
    xlabel('Frequency [Hz]')
    ylabel('Accelerance [m/s$^2$/N]')
%     title(casestoread{t})
end

%% EXPORT RESULTS TO .CSV FILE

for i=1:length(casestoread)
    T(i,:) = table(casestoread(i),beam.m(i),beam.L(i),beam.b(i),beam.t(i),beam.rho(i),test.(casestoread{i}).rf',test.(casestoread{i}).af', ...
        test.(casestoread{i}).E_iso',test.(casestoread{i}).E_af',test.(casestoread{i}).E_rf');
end
T = renamevars(T,["Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11"], ...
    ["Beam","Mass (kg)","Length (m)","Width (m)","Thickness (m)","Desnity (kg/m3)","Resonance Frequencies (Hz)","Antirresonance Frequencies (Hz)","E_ISO (GPa)","E_AF (GPa)","E_RF (GPa)"]);

% filename = 'Set_test3.xlsx';              % Rename the file name when changing the data files
% writetable(T,filename,"WriteMode","append","PreserveFormat",true)
%% FUNCTIONS

function E_I = E_ISO(L,rhom,Ixx,f)                           % Young's modulus with ISO-16940 method
    lambda = [1.87510 4.69410 7.85476];
    E_I = rhom/Ixx*(2*pi*(L/2)^2*f./lambda.^2).^2*1E-9;
end

function E_A = E_AFF(L,rhom,Ixx,f)                           % Young's modulus with antiresonance frequencies of a free-free beam method
    lambda = [3.75038 9.39740 15.73438];
    E_A = rhom/Ixx*(2*pi*L^2*f./lambda.^2).^2*1E-9;
end

function E_R = E_RFF(L,rhom,Ixx,f)                           % Young's modulus with resonance frequencies of a free-free beam method
    lambda = [4.73004 10.99561 17.27876];
    E_R = rhom/Ixx*(2*pi*L^2*f./lambda.^2).^2*1E-9;
end
