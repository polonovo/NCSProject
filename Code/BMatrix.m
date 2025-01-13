% Construct Beta matrix using Japanese data

clear, close all

%% Import data

addpath("./../Data")
data = readtable("japanData2.csv");

% get rid of Kumamoto and extra jKpan
adj = data(3:49,5:end);
adj(:,"Kumamoto") = [];
adj(45,:) = [];
adj = table2array(adj)'./100;

adj(ad)

% Make row stochastic
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

% Save prefecture names
pref_names = data.Var2(3:49);
pref_names(43) = [];

n = length(adj);

save("Beta.mat",'adj','pref_names')
