clc
clear
%% load FC-matrix data
path1 = 'F:\AAL\1';
path2 = 'F:\AAL\2';
path3 = 'F:\AAL\3';

list1 = dir([path1,filesep,'*.mat']);
list2 = dir([path2,filesep,'*.mat']);
list3 = dir([path3,filesep,'*.mat']);

z1 = [];
z2 = [];
z3 = [];

for i = 1:length(list1)
    file1 = [path1,filesep,list1(i).name]
    file_load = load(file1);
    file_mat = file_load.ROICorrelation_FisherZ;
    z1(:,:,i) = file_mat;
end
for i = 1:length(list2)
    file2 = [path2,filesep,list2(i).name]
    file_load = load(file2);
    file_mat = file_load.ROICorrelation_FisherZ;
    z2(:,:,i) = file_mat;
end
for i = 1:length(list3)
    file3 = [path3,filesep,list3(i).name]
    file_load = load(file3);
    file_mat = file_load.ROICorrelation_FisherZ;
    z3(:,:,i) = file_mat;
end
%% exclude infinite value
z1(z1 == Inf) = 0;
z2(z2 == Inf) = 0;
z3(z3 == Inf) = 0;
%% calculate SEVERE-by-MILD-baseline and SEVERE-by-MILD-4wk correlation matrix
AA1M = cat(3,z2,z3);
[x1,y1,n1] = size(z1);
[x2,y2,n2] = size(AA1M);
acc = [];
for i = 1:n1
    subj = z1(:,:,i);
    for j = 1:n2
        subj_c = AA1M(:,:,j);
        cor = corrcoef(subj,subj_c);
        rvalue = cor(1,2);
        acc(i,j) =  rvalue;
    end
end
%% Calculate the ratio of subjects who fitted our hypothesis
SMb = acc(:,1:21);
SM4wk = acc(:,22:end);
final_acc = [];
final_acc_mean=[];
for i = 1:n1
    SMb_subj = SMb(i,:);
    SM4wk_subj = SM4wk(i,:);
     meanvalue1 = mean(SMb_subj);
     meanvalue2 = mean(SM4wk_subj);
  
 if meanvalue1>=meanvalue2
        final_acc(i) = 0
    else 
            final_acc(i) = 1
    end
 
end
final_ratio = mean(final_acc);


%% Permutation test
no_iteration = 5000; % Permutation times
prediction_r = zeros(no_iteration,1);
prediction_r(1,1) = final_ratio; % Real ratio

for it=2:no_iteration
    
    fprintf('\n Performing iteration %d out of %d',it, no_iteration);


final_acc = [];
for i = 1:n1
    radp = randperm(39);
    SEVERE_to_all = acc(i,:);
    SEVERE_to_all_rand = SEVERE_to_all(radp);
    SMb_subj = SEVERE_to_all_rand(1:21);
    SM4wk_subj = SEVERE_to_all_rand(22:end);
    meanvalue1 = mean(SMb_subj);
    meanvalue2 = mean(SM4wk_subj);

    if  meanvalue1 >=  meanvalue2
        final_acc_mean(i) = 0;
    else 
            final_acc_mean(i) = 1;
    end
end

final_ratio = mean(final_acc_mean);
prediction_r(it) = final_ratio;
end

%% Calculate the P-value of permutation test
pval_r = mean((prediction_r(2:end) > prediction_r(1)));



