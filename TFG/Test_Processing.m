% ***********************************************************************
%                   DATA PROCESSING OF DYNAMIC TESTING
% ***********************************************************************

% Collects all the numeric data from the plate's testing and places them
% into different structures to automate the post-processing

clc;clear;close all

%% BEAM'S PROPERTIES

plate.input = importdata('Ensayos_Placa_Ref\geomdata.txt',';',2);               % CHANGE THE FILE
plate.m = plate.input.data(:,1);
plate.L = plate.input.data(:,2);
plate.b = plate.input.data(:,3);
plate.t = plate.input.data(:,4);
plate.Ixx = plate.b.*plate.t.^3./12;                      
plate.rho = plate.m./plate.L./plate.b./plate.t;             
plate.rhom = plate.rho.*plate.b.*plate.t;                 

%% CREATE TREE STRUCTURE

% 3-D printing parameters

plate.seed{1} = {'F','V','E'};                                       % Printing orientation (Flat, Vertical, On Edge)
plate.seed{2} = {'99','90' '85','70','60','50','40'};                % Infill percentage (%) 
plate.seed{3} = {'00','15','30','45','60','90'};                     % Raster angle (deg)
plate.seed{4} = {'1','2','3','4'};                                   % Layer thickness (1E-4 m)
plate.seed{5} = {'F','S'};                                           % Number of test (Firts, Second)   
plate.seed{6} = {'a','b','c','d'};                                   % Sets of the same beam (4 test of each one)

 % Specify the files that are going to be read here
casestoread = {'F99001Fa', 'F99001Fb', 'F99001Fc', 'F99001Fd', 'F99001Sa', 'F99001Sb', 'F99001Sc', 'F99001Sd',
};


% Generate all possible combinations of beams

cont = 1;
for i=1:length(plate.seed{1})
    for j=1:length(plate.seed{2})
        for k=1:length(plate.seed{3})
            for l=1:length(plate.seed{4})
                for m=1:length(plate.seed{5})
                    for n=1:length(plate.seed{6})
                        plate.name{cont,1} = [plate.seed{1}{i} plate.seed{2}{j} plate.seed{3}{k} plate.seed{4}{l} plate.seed{5}{m} plate.seed{6}{n}];
                        cont = cont + 1;
                    end
                end
            end
        end
    end
end

for s=1:length(casestoread)
    filecontent = readmatrix(['Ensayos_Placa_Ref\' casestoread{s} '.csv'],'Delimiter',';','NumHeaderLines',20,'DecimalSeparator',',');             % CHANGE THE FILE
    test.(casestoread{s}).freq = filecontent(:,1);
    test.(casestoread{s}).real = filecontent(:,2);
    test.(casestoread{s}).imag = filecontent(:,3);
    
    test.(casestoread{s}).amp = sqrt(filecontent(:,2).^2 + filecontent(:,3).^2);

end

%% FIGURES
figure()
tiledlayout(4,2)
for t=1:length(casestoread)
    nexttile
    hold on
    plot(test.(casestoread{t}).freq,log(test.(casestoread{t}).amp))
    title(casestoread{t})
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