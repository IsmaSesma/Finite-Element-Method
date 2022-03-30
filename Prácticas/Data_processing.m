% *************************************************************
%            DATA PROCESSING OF DYNAMIC TESTING
% *************************************************************

% Collects all the numeric data from the beam's testing and places them
% into different structures to automate the post-processing

clc;clear;close all

s = 50;                                                 % Number of cases to analyze

for i = 1:s
cases.names{i} = ['caso' num2str(i) '.txt'];
end



for i = 1:s
    ensayos.(cases.names{i}).fr=i;
end


for caso = 10:20
    ensayos.(cases.names{caso}).fr2=ensayos.(cases.names{caso}).fr^2;
end