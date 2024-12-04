% Construct Beta matrix using Japanese data

clear, close all

%% Import data

addpath("./../Data")
data = readtable("japanData2.csv");

adj = data(3:49,5:end);
adj(:,"Kumamoto") = [];
adj(45,:) = [];
adj = table2array(adj)'./100;

deg = sum(adj,2);
adj = (diag(1./deg))*adj;

% Make matrix double stochastic

% iter = 40;
% 
% for i = 1:iter
%     deg = sum(adj,2);
%     adj = (diag(1./deg))*adj;
%     deg2 = sum(adj,1);
%     adj = adj*(diag(1./deg2));
% end

adj = 0.5*adj;

names = data.Var2(3:49);
names(43) = [];

n = length(adj);

save("Beta.mat",'adj')
