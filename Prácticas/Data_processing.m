% ***********************************************************************
%                   DATA PROCESSING OF DYNAMIC TESTING
% ***********************************************************************

% Collects all the numeric data from the beam's testing and places them
% into different structures to automate the post-processing

clc;clear;close all

%% BEAM'S PROPERTIES

beam.input = importdata('data_test1\geomdata.txt',';',2);               % CHANGE THE FILE
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
casestoread = {'E99001ZZa',
'E90451ZZa',
'E85001ZZa',
'E70001Ta',
'E60001Ta',
'E50001Ta',
'E40001Ta',
'E99002ZZa',
'E99003ZZa',
'E99004ZZa',
'E99001ZZb',
'E90451ZZb',
'E85001ZZb',
'E70001Tb',
'E60001Tb',
'E50001Tb',
'E40001Tc',
'E99002ZZb',
'E99003ZZb',
'E99004ZZb',
'E99001ZZc',
'E90451ZZc',
'E85001ZZc',
'E70001Tc',
'E60001Tc',
'E50001Tc',
'E40001Tc',
'E99002ZZc',
'E99003ZZc',
'E99004ZZc'
};


% Generate all possible combinations of beams

cont = 1;
for i=1:length(beam.seed{1})
    for j=1:length(beam.seed{2})
        for k=1:length(beam.seed{3})
            for l=1:length(beam.seed{4})
                for m=1:length(beam.seed{5})
                    for n=1:length(beam.seed{6})
                        beam.name{cont,1} = [beam.seed{1}{i} beam.seed{2}{j} beam.seed{3}{k} beam.seed{4}{l} beam.seed{5}{m} beam.seed{6}{n}];
                        cont = cont + 1;
                    end
                end
            end
        end
    end
end

for s=1:length(casestoread)
    filecontent = readmatrix(['data_test1\' casestoread{s} '.csv'],'Delimiter',';','NumHeaderLines',20,'DecimalSeparator',',');             % CHANGE THE FILE
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
tiledlayout(5,6)
for t=1:length(casestoread)
    nexttile
    hold on
    plot(test.(casestoread{t}).freq,log(test.(casestoread{t}).amp))
    plot(test.(casestoread{t}).rf,test.(casestoread{t}).peak,'ro')
    plot(test.(casestoread{t}).af,test.(casestoread{t}).min,'bo')
    title(casestoread{t})
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
