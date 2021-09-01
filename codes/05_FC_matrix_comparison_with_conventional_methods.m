clc
clear
%% Load FC-matrix data
path1 = 'D:\Fragmented_causal_V2\4-GSR\AAL_globalCF\C';
path2 = 'D:\Fragmented_causal_V2\4-GSR\AAL_globalCF\A';
path3 = 'D:\Fragmented_causal_V2\4-GSR\AAL_globalCF\A1M';
list1 = dir([path1,filesep,'*.mat']);
list2 = dir([path2,filesep,'*.mat']);
list3 = dir([path3,filesep,'*.mat']);
% Only keep the upper triangle
NumROI = 116;
conn_msk = ones(116);
Ind_01    = find(triu(ones(116),1));
Ind_02    = find(conn_msk(Ind_01) ~= 0);

for i = 1:length(list1)
    file1 = [path1,filesep,list1(i).name]
    file_load = load(file1);
    file_mat = file_load.ROICorrelation_FisherZ;
    DATA_SEVERE(i,:) = file_mat(Ind_01(Ind_02));
end
for i = 1:length(list2)
    file2 = [path2,filesep,list2(i).name]
    file_load = load(file2);
    file_mat = file_load.ROICorrelation_FisherZ;
    DATA_MILD_baseline(i,:) = file_mat(Ind_01(Ind_02));
end
for i = 1:length(list3)
    file3 = [path3,filesep,list3(i).name]
    file_load = load(file3);
    file_mat = file_load.ROICorrelation_FisherZ;
    DATA_MILD_4wk(i,:) = file_mat(Ind_01(Ind_02));
end

%% 2-sample t-test
[~,p1,~,stat1] = ttest2(DATA_SEVERE,DATA_MILD_baseline);
mask = p1<0.001; %build a mask
t1 = stat1.tstat;
t1 = t1.*mask;
p1 = p1.*mask;

[~,p2,~,stat2] = ttest2(DATA_SEVERE,DATA_MILD_4wk);
t2 = stat2.tstat;
t2 = t2.*mask;
p2 = p2.*mask;

[~,p3,~,stat3] = ttest2(DATA_MILD_4wk,DATA_MILD_baseline);
t3 = stat3.tstat;
t3 = t3.*mask;
p3 = p3.*mask;

%% Mean changes of selected FCs after 4-week social media use
m1 = (mean(DATA_MILD_4wk,1)-mean(DATA_MILD_baseline,1)).*mask;
fcm = zeros(116);
fcm(Ind_01(Ind_02))=m1;

%% Exclude FCs outside the mask
DATA_ALL = cat(1,DATA_SEVERE,DATA_MILD_baseline,DATA_MILD_4wk);
DATA_ALL2 = DATA_ALL.*repmat(mask,68,1);

j = 1
for i = 1:6670
    lie = DATA_ALL2(:,i);
    
    if sum(lie) ==0
        kk(j) = i;
        j = j+1;
    else
       
    end
end
DATA_ALL2(:,kk) = [];

%% The comparison of selected FCs between 3 groups
SEVERE_group = DATA_ALL2(1:29,:);
MILD_baseline_group = DATA_ALL2(30:50,:);
MILD_4wk_group = DATA_ALL2(51:end,:);