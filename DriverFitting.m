% Optimization Driver/ Sensitivity analysis
% clear; 
%    close all; clc;
tic
% plotting while the optimization runs
flag_plot = 1; % flag_plot = 1 ==>  plot
%%%%%%%Producing the first Generetion%%%%%%%%
%% adjustable variables are multiplied to the old variables

%var ={ 1.64860037641656	1.30776963814884	1.45838224315620
%0.766604755794347	1.58780058107716	1.66979993127668	1.43361221916716	1.63061398905641	1.24984502631487	1.95656814763700	1.92727542609161	1.01532003606808	0.716314335317795	0.672539492283114	1.65502459182752	1.78304757052419]
% var = [   1.8616/1.5    1.7228    1.9432    0.9809    1.0248    0.8868    1.4304*1.5    1.5103    0.8826    4.4748*1.1    2.1132    0.9160    1.0171    0.9865    1.9299    1.0519]
var = [   1.8616/1.5    1.7228    1.9432  0.9809    1.0248    0.8868    1.4304*1.5    1.5103    0.8826    4.4748*1.1    2.1132    1    1    1    1.9299    1*0.1]

[hill_19 Ca50_19 hill_23 Ca50_23]= getpCa50 (var);       

Max_Force= getMaxForce (var) 


Data.freq =[2 3 4 5 6 7 8 9 10 4 4 4 ];
filename = {'Ca_2Hz.mat', 'Ca_3Hz.mat', 'Ca_4Hz.mat', 'Ca_5Hz.mat', 'Ca_6Hz.mat', 'Ca_7Hz.mat', 'Ca_8Hz.mat', 'Ca_9Hz.mat', 'Ca_10Hz.mat','Ca_4Hz.mat','Ca_4Hz.mat','Ca_4Hz.mat'};
% Set temperature fot the experiment environment
Data.TmpC = 37; % centigrade 
err  = objective_fun_rat_data(var,Data,flag_plot)
toc