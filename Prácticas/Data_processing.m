% ***********************************************************************
%                   DATA PROCESSING OF DYNAMIC TESTING
% ***********************************************************************

% Collects all the numeric data from the beam's testing and places them
% into different structures to automate the post-processing

clc;clear;close all

%% BEAM'S PROPERTIES

beam.input = importdata('data\geomdata.txt',';',2);
beam.m = beam.input.data(:,1);
beam.L = beam.input.data(:,2);
beam.b = beam.input.data(:,3);
beam.t = beam.input.data(:,4);
beam.Ixx = beam.b.*beam.t.^3./12;                      
beam.rho = beam.m./beam.L./beam.b./beam.t;             
beam.rhom = beam.rho.*beam.b.*beam.t;                 

%% CREATE TREE STRUCTURE

% 3-D printing parameters

test.seed{1} = {'F','V','E'};                                       % Printing orientation (Flat, Vertical, On Edge)
test.seed{2} = {'99','90' '85','70','60','50','40'};                % Infill percentage (%) 
test.seed{3} = {'00','15','30','45','60','90'};                     % Raster angle (deg)
test.seed{4} = {'1','2','3','4'};                                   % Layer thickness (1E-4 m)
test.seed{5} = {'ZZ','T'};                                          % Infill pattern (Zig-Zag, Triangular)   
test.seed{6} = {'a','b','c'};                                       % Sets of the same beam (3 models of each one)

casestoread = {'E99001ZZa','E90451ZZa','E85001ZZa','E70001Ta','E60001Ta','E50001Ta','E40001Ta','E99002ZZa','E99003ZZa','E99004ZZa',...
               'E99001ZZb','E90451ZZb','E85001ZZb','E70001Tb','E60001Tb','E50001Tb','E99002ZZb','E99003ZZb','E99004ZZb',...
               'E99001ZZc','E90451ZZc','E85001ZZc','E70001Tc','E60001Tc','E50001Tc','E99002ZZc','E99003ZZc','E99004ZZc'};    % Specify the files that are going to be read here

% Generate all possible combinations of beams

cont = 1;
for i=1:length(test.seed{1})
    for j=1:length(test.seed{2})
        for k=1:length(test.seed{3})
            for l=1:length(test.seed{4})
                for m=1:length(test.seed{5})
                    for n=1:length(test.seed{6})
                        test.name{cont,1} = [test.seed{1}{i} test.seed{2}{j} test.seed{3}{k} test.seed{4}{l} test.seed{5}{m} test.seed{6}{n}];
                        cont = cont + 1;
                    end
                end
            end
        end
    end
end

for s=1:length(casestoread)
    filecontent = readmatrix(['data\' casestoread{s} '.csv'],'Delimiter',';','NumHeaderLines',20,'DecimalSeparator',',');
    test.(casestoread{s}).freq = filecontent(:,1);
    test.(casestoread{s}).real = filecontent(:,2);
    test.(casestoread{s}).imag = filecontent(:,3);
    
    test.(casestoread{s}).amp = sqrt(filecontent(:,2).^2 + filecontent(:,3).^2);
    
    % MCM 2022-04-25
    [test.(casestoread{s}).peak, test.(casestoread{s}).rf] = findpeaks(test.(casestoread{s}).amp,test.(casestoread{s}).freq,'MinPeakProminence',50,'NPeaks',3);          % Find the resonance frequencies
    tolook_f = test.(casestoread{s}).freq(10:end);
    tolook_v = 1./test.(casestoread{s}).amp(10:end);
    [test.(casestoread{s}).min, test.(casestoread{s}).af] = findpeaks(tolook_v,tolook_f,'MinPeakProminence',0.065,'NPeaks',3);          % Find the antiresonance frequencies
    test.(casestoread{s}).min = 1./test.(casestoread{s}).min;
    % ----

    test.(casestoread{s}).E_iso = (E_ISO(beam.L(s),beam.rhom(s),beam.Ixx(s),test.(casestoread{s}).af));
    test.(casestoread{s}).E_af = diag(E_AFF(beam.L(s),beam.rhom(s),beam.Ixx(s),test.(casestoread{s}).af));
    test.(casestoread{s}).E_rf = diag(E_RFF(beam.L(s),beam.rhom(s),beam.Ixx(s),test.(casestoread{s}).rf));

end

%%
% MCM 2022-04-25
figure()
tiledlayout(6,4)
for t=1:28
    nexttile
    hold on
    plot(test.(casestoread{t}).freq,test.(casestoread{t}).amp)
    plot(test.(casestoread{t}).rf,test.(casestoread{t}).peak,'ro')
    plot(test.(casestoread{t}).af,test.(casestoread{t}).min,'bo')

end

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