% Creat delta matrix from hospital data

clear,close all

addpath('./../Data')

%% Import data

load('hospitalsper100k.mat')
migata = readtable('japanData2.csv');

data = a109I;

% Eliminate some rows (Kumamoto, mean, std)
data([1 44 48 49],:) = [];

data = table2array(data(:,2));

mu = mean(data); % Mean

datan = data./mu; % Normalize with mean so valuec lose to average = 1

% Exponential function to assign recovery rate based on number of hospitals
% per 100k inhabitants. Try to make average = 0.5 recovery rate
% f = @(x) exp(x/2.71-1)-0.2;
f = @(x) exp(x/2.71-1.2);

Delta = f(datan);

plot(Delta)

save('Delta.mat','Delta')

