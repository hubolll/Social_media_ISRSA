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

%% Load FC matrix
label = load('FC_matrix');
    
%% Convert FC matrix to 3-conlumns matrix
label2 = []; k=1;
for i=1:116
    for j=1:116
        if label(i,j) ~= 0
            label2(k,1) = net(i,1);
            label2(k,2) = net(j,1);
            label2(k,3) = label(i,j);
             k=k+1;
        else
            label2 = label2;
        end
    end
end

%% Convert FCs between brain regions to FCs between brain networks
Net_positive = zeros(9,9);
Net_negative = zeros(9,9);
p = 0
n = 0
[m,n] = size(label2);
for i = 1:m
    x = label2(i,1)
    y = label2(i,2)
    z = label2(i,3)
    label_a = find(net(:,1)==x);
    neta = net(label_a,2);
    lable_b = find(net(:,1)==y);
    netb = net(lable_b,2);
    % Separate positive and negative FCs
    if z>0
        Net_positive(neta,netb) = Net_positive(neta,netb)+z
    else 
        Net_negative(neta,netb) = Net_negative(neta,netb)+z
    end
end

%% Convert FC matrix between brain networks to 3-conlumns matrix
NumROI = 9;
Net_positive_4 = zeros(9,9);
conn_msk = ones(9);
Ind_01    = find(triu(ones(9),0));
Ind_02    = find(conn_msk(Ind_01) ~= 0);
Net_positive_3 = Net_positive(Ind_01(Ind_02));
Net_positive_4(Ind_01(Ind_02)) = Net_positive_3;
label3 = []; k=1;
for i=1:9
    for j=1:9
        if Net_positive_4(i,j) ~= 0
            label3(k,:) = [i j Net_positive_4(i,j)]; k=k+1;
        else
            label3 = label3;
        end
    end
end


