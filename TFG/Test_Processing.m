% ***********************************************************************
%                   DATA PROCESSING OF DYNAMIC TESTING
% % **********************************************************************

% Collects all the numeric data from the plate's testing and places them
% into different structures to automate the post-processing

clc;clear;close all

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth',2)

f = 1000; vf = (1:f);

% Specify the files that are going to be read here
casestoread1 = {'Placa_Ref', 'Placa_250_acel_bot'};
casestoread2 = {'Bot_Ref_Data_FEM','Bot_Holed_Data_FEM', 'Bot_180Hz_Data_FEM', 'Bot_250Hz_Data_FEM'};

for s=1:length(casestoread1)
    filecontent = readmatrix(['Plate_Tunned_Test/' casestoread1{s} '.csv'],'Delimiter',';','NumHeaderLines',20,'DecimalSeparator',',');          
    test.(casestoread1{s}).freq = filecontent(:,1);
    test.(casestoread1{s}).real = filecontent(:,2);
    test.(casestoread1{s}).imag = filecontent(:,3);
    
    test.(casestoread1{s}).amp = sqrt(filecontent(:,2).^2 + filecontent(:,3).^2);
end

for s=1:length(casestoread2)
    filecontent = readmatrix(['Plate_Tunned_Test/' casestoread2{s} '.csv'],'DecimalSeparator','.');            
    
    test.(casestoread2{s}).amp = abs(filecontent(:,1));
end

%% FIGURES

figure()
tiledlayout(1,1)
nexttile(1)
hold on
%      plot(test.(casestoread1{1}).freq,test.(casestoread1{1}).amp)
%      plot(test.(casestoread1{2}).freq,test.(casestoread1{2}).amp)
%     plot(test.(casestoread1{3}).freq,test.(casestoread1{3}).amp)
%     semilogy(vf,test.(casestoread2{1}).amp)
    semilogy(vf,test.(casestoread2{2}).amp)
    semilogy(vf,test.(casestoread2{3}).amp)
%     semilogy(vf,test.(casestoread2{4}).amp)
    set(gca,'XScale','lin','YScale','log')
    set(gca,'XLim',[0,400])
    set(gca, 'YLim', [10^-3, 10^4])
    box;grid;
    legend('Holed plate','Tuned plate, model II' ,'Location','SouthEast')
    xlabel('Frequency (Hz)')
    ylabel('Accelerance [m/s$^2$/N]')

% Comparison between reference and tunned plates
% figure()
% plot(test.Placa_1.freq,(test.Placa_tunned_3_2.amp))
% hold on
% plot(test.Placa_Ref.freq,(test.Placa_Ref.amp))
% set(gca,'YScale','log')
% set(gca,'XLim',[0, 400])
% xlabel('Frequency [Hz]')
% ylabel('Amplitude [m/s$^2$]')
% legend('Tunned Plate','Reference Plate')
