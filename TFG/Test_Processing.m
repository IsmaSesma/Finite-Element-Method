% ***********************************************************************
%                   DATA PROCESSING OF DYNAMIC TESTING
% % **********************************************************************

% Collects all the numeric data from the plate's testing and places them
% into different structures to automate the post-processing

clc;clear;close all

% Specify the files that are going to be read here
casestoread = {'Ensayo_Pletina_Blanca', 'Placa_1_ac1', 'Placa_1', 'Placa_2_ac1', 'Placa_1_ac2', 'Placa_Ref','Placa_tunned_3_2'};

for s=1:length(casestoread)
    filecontent = readmatrix(['Plate_Tunned_Test/' casestoread{s} '.csv'],'Delimiter',';','NumHeaderLines',20,'DecimalSeparator',',');             % CHANGE THE FILE
    test.(casestoread{s}).freq = filecontent(:,1);
    test.(casestoread{s}).real = filecontent(:,2);
    test.(casestoread{s}).imag = filecontent(:,3);
    
    test.(casestoread{s}).amp = sqrt(filecontent(:,2).^2 + filecontent(:,3).^2);

end

%% FIGURES
figure()
tiledlayout(3,3)
for t=1:length(casestoread)
    nexttile
    hold on
    plot(test.(casestoread{t}).freq,log(test.(casestoread{t}).amp))
    title(casestoread{t})
end

% Comparison between reference and tunned plates
figure()
plot(test.Placa_1.freq,(test.Placa_tunned_3_2.amp))
hold on
plot(test.Placa_Ref.freq,(test.Placa_Ref.amp))
set(gca,'YScale','log')
set(gca,'XLim',[0, 1000])
xlabel('Frequency [Hz]')
ylabel('Amplitude [m/s$^2$]')
legend('Tunned Plate','Reference Plate')
