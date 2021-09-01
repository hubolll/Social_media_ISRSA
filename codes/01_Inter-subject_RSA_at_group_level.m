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
%% Comparing the SEVERE-MILD-baseline similarity to SEVERE-MILD-4wk similarity
SEVERE-MILD-baseline = acc(:,1:21);
SEVERE-MILD-4wk = acc(:,22:end);
SMb = SEVERE-MILD-baseline(:);
SM4wk = SEVERE-MILD-4wk(:);
[~,p,~,stat] = ttest2(SMb,SM4wk);
mean(SMb)
mean(SM4wk)
std(SMb)
std(SM4wk)
