% ***********************************************************************
%                   DATA PROCESSING OF DYNAMIC TESTING
% ***********************************************************************

% Collects all the numeric data from the beam's testing and places them
% into different structures to automate the post-processing

clc;clear;close all

%% BEAM'S PROPERTIES

beam.m = 0.03;
beam.L = 0.2;
beam.b = 0.02;
beam.t = 0.004;
beam.Ixx = beam.b*beam.t^3/12;                      
beam.rho = beam.m/beam.L/beam.b/beam.t;             
beam.rhom = beam.rho*beam.b*beam.t;                 

%% CREATE TREE STRUCTURE

% 3-D printing parameters

test.seed{1} = {'F','V','E'};                                       % Printing orientation (Flat, Vertical, On Edge)
test.seed{2} = {'99','90' '85','70','60','50','40'};                % Infill percentage (%) 
test.seed{3} = {'00','15','30','45','60','90'};                     % Raster angle (deg)
test.seed{4} = {'1','2','3','4'};                                   % Layer thickness (1E-4 m)
test.seed{5} = {'ZZ','T'};                                          % Infill pattern (Zig-Zag, Triangular)   
test.seed{6} = {'a','b','c'};                                       % Sets of the same beam (3 models of each one)

casestoread = {'F90001ZZa'};                  % Specify the files that are going to be read here

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
    filecontent = importdata([casestoread{s} '.txt'],',');
    test.(casestoread{s}).freq = filecontent(:,1);
    test.(casestoread{s}).amp = filecontent(:,2);
    test.(casestoread{s}).headers = filecontent;

    [test.(casestoread{s}).peak, test.(casestoread{s}).rf] = findpeaks(test.(casestoread{s}).amp);          % Find the resonance frequencies
    [test.(casestoread{s}).min, test.(casestoread{s}).af] = findpeaks(-test.(casestoread{s}).amp);          % Find the antiresonance frequencies

    test.(casestoread{s}).E_iso = E_ISO(beam.L,beam.rhom,beam.Ixx,test.(casestoread{s}).af);
    test.(casestoread{s}).E_af = E_AFF(beam.L,beam.rhom,beam.Ixx,test.(casestoread{s}).af);
    test.(casestoread{s}).E_rf = E_RFF(beam.L,beam.rhom,beam.Ixx,test.(casestoread{s}).rf);

end


% for s=1:2
%     filecontent = importdata([test.name{s} '.txt'],'\t',3);
%     test.(test.name{s}).t = filecontent.data(:,1);
%     test.(test.name{s}).Amp = filecontent.data(:,2);
%     test.(test.name{s}).names = filecontent.colheaders;
% end



% figure('Color','White')
% hold on
% leg = '';
% for s=1:2
%     plot(test.(test.name{s}).t,test.(test.name{s}).Amp)
%     leg{s} = test.name{s};
% end

%% FUNCTIONS

function E_I = E_ISO(L,rhom,Ixx,f)                           % Young's modulus with ISO-6940 method
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