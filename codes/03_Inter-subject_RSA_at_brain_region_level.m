clc
clear

%% load FC-matrix data
path1 = 'F:\Fragmented_New_test\6-cause_NEW_test_rFC\AAL\C';
path2 = 'F:\Fragmented_New_test\6-cause_NEW_test_rFC\AAL\A';
path3 = 'F:\Fragmented_New_test\6-cause_NEW_test_rFC\AAL\A1M';

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
%% calculate SEVERE-by-MILD-baseline and SEVERE-by-MILD-4wk correlation matrix of each brain region
for all_it = 1:116 %116 brain regions in total
    FC1 = z1(all_it,:,:);
    FC2 = z2(all_it,:,:);
    FC3 = z3(all_it,:,:);
FC_AA1M = cat(3,FC2,FC3);
[x1,y1,n1] = size(FC1);
[x2,y2,n2] = size(FC_AA1M);


for i = 1:n1
    subj = FC1(:,:,i);
    for j = 1:n2
        subj_c = FC_AA1M(:,:,j);
        cor = corrcoef(subj,subj_c);
        rvalue = cor(1,2);
        acc(i,j,all_it) =  rvalue;
    end
end
%% Calculate the ratio of subjects who fitted our hypothesis
SMb = acc(:,1:21,all_it);
SM4wk = acc(:,22:end,all_it);
final_acc = [];
final_acc_mean=[];
for i = 1:n1
    SMb_subj = SMb(i,:);
    SM4wk_subj = SM4wk(i,:);
    meanvalue1 = mean(SMb_subj);
    meanvalue2 = mean(SM4wk_subj);

    if meanvalue1>meanvalue2
        final_acc(i) = 0;
        final_acc_verse(i) = 1;
    else 
            final_acc(i) = 1;
            final_acc_verse(i) = 0;
    end
  
end
final_ratio = mean(final_acc); % Ratio of subjects who fitted our hypothesis
final_ratio_verse = mean(final_acc_verse); % Ratio of subjects who was opposed to our hypothesis
ROI_ratio(all_it) = final_ratio;
ROI_ratio_verse(all_it) = final_ratio_verse;
end

%% Permutation test
no_iteration = 5000;
for all_it =1:116 % Permutation test for each brain region
prediction_r(:,all_it) = zeros(no_iteration,1);
prediction_r2(:,all_it) = zeros(no_iteration,1);

prediction_r(1,all_it) = ROI_ratio(all_it);
prediction_r2(1,all_it) = ROI_ratio_verse(all_it);
for it=2:no_iteration

fprintf('\n Performing iteration %d out of %d',it, no_iteration);   

for i = 1:n1
    
radp = randperm(39);
    rand_acc = acc(i,radp,all_it);
  SMb_subj = rand_acc(1:21);
    SM4wk_subj = rand_acc(22:end);
    meanvalue1 = mean(SMb_subj);
    meanvalue2 = mean(SM4wk_subj);

    if  meanvalue1 >  meanvalue2
        final_acc_mean(i) = 0;
         final_acc_mean_verse(i) = 1;
    else 
            final_acc_mean(i) = 1;
            final_acc_mean_verse(i) = 0;
    end

end
final_ratio = mean(final_acc_mean);
final_ratio_verse = mean(final_acc_mean_verse);

prediction_r(it,all_it) = final_ratio;
prediction_r2(it,all_it) = final_ratio_verse;

end

%% Calculate the P-value of permutation test

pval_r(all_it) = mean((prediction_r(2:end,all_it) > prediction_r(1,all_it)));
pval_r2(all_it) = mean((prediction_r2(2:end,all_it) > prediction_r2(1,all_it)));
end



