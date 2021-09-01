clc
clear
%% Load FC_network_conversion file
net = load('D:\Fragmented_article\ÕýÊ½\20201026_NBS\net_conversion.txt');
% This file could convert AAL altas labels into Schaefer's 7-Networks atlas
% (2 additional networks are basal ganglia and cerebellum)
% 1-Visual network, 2-Somato-motor network, 3-Dorsal attention network
% 4-Ventral attention network, 5-Limbic system, 6-Fronto-parietal network
% 7-Default mode network, 8-Basal ganglia, 9-Cerebellum
% We used a winner-take-all strategy to make this conversion file

%% Load brain regions
    label = ; % P-value of selected brain regions

%% Assign brain regions with significance into brain networks
    k = [];
for i = 1:116
    region = label(i);
    if region<0.05
      network = net(i,2);
      k(network) = k(network)+1;
    else
        network = net(i,2);
        k(network) = k(network)
    end
end

