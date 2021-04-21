%% NR-F enhances the quality of Hi-C contact map from the human genome
%  This is an example on chromosome 19


%%
clear
clc
addpath(genpath(pwd));

%% load data
% global parameter
resolution = '100kb';
downsample_ratio = 1/100;
m = 2;
error1 = zeros(5,1);
error2 = zeros(5,1);
error3 = zeros(5,1);
loops = 5;
cutoff = 0.95;
tic
%%
for loop = 1:loops
    error1_old = error1;
    error2_old = error2;
    error3_old = error3;
    
    for chromosome = 19
        
        [contact_mat,downsamped_contact_mat] = Import_paired_Data(['chr',num2str(chromosome)],...
            resolution,downsample_ratio);%paired contact mat
        cluster = tabulate(louvain(contact_mat));
        sub_window_size = round(max(cluster(:,2)));
        
        %% ***************************************************************
        
        disp('NR-F Denoising.....')
        NR_F_contact_mat = denoise(downsamped_contact_mat,sub_window_size,m,1);
        
        disp('NR-F twice Denoising.....')
        NR_F_contact_mat_twice = denoise(NR_F_contact_mat,sub_window_size,m,1);
        
        disp('NE Denoising.....')
        NE_contact_mat = denoise(downsamped_contact_mat,sub_window_size,m,2);
        disp('ND Denoising.....')
        ND_contact_mat = denoise(downsamped_contact_mat,sub_window_size,m,3);
        
      
        %% ***********************************************************************
        cutoff1 = quantile(contact_mat(:),cutoff);
        contact_mat(contact_mat>=cutoff1) = cutoff1;
        
        cutoff2 = quantile(downsamped_contact_mat(:),cutoff);
        downsamped_contact_mat(downsamped_contact_mat>=cutoff2) = cutoff2;
        
        cutoff3 = quantile(NR_F_contact_mat(:),cutoff);
        NR_F_contact_mat(NR_F_contact_mat>=cutoff3) = cutoff3;
        
        cutoff4 = quantile(NE_contact_mat(:),cutoff);
        NE_contact_mat(NE_contact_mat>=cutoff4) = cutoff4;
        
        cutoff5 = quantile(ND_contact_mat(:),cutoff);
        ND_contact_mat(ND_contact_mat>=cutoff5) = cutoff5;
        
        cutoff6 = quantile(NR_F_contact_mat_twice(:),cutoff);
        NR_F_contact_mat_twice(NR_F_contact_mat_twice>=cutoff6)= cutoff6;
        
        %% **********************************************************************
        downsamped_contact_mat = downsamped_contact_mat/max(max(downsamped_contact_mat));
        contact_mat = contact_mat/max(max(contact_mat));
        NR_F_contact_mat = NR_F_contact_mat/max(max(NR_F_contact_mat));
        NR_F_contact_mat_twice = NR_F_contact_mat_twice/max(max(NR_F_contact_mat_twice));
        NE_contact_mat = NE_contact_mat/max(max(NE_contact_mat));
        ND_contact_mat = ND_contact_mat/max(max(ND_contact_mat));
        
        
        %% ***********************************************************************
        disp('Evaluating.....')
        error1(1) = ssim_sub(downsamped_contact_mat,contact_mat,sub_window_size);
        error1(2) = ssim_sub(NR_F_contact_mat,contact_mat,sub_window_size);
        error1(3) = ssim_sub(NR_F_contact_mat_twice,contact_mat,sub_window_size);
        error1(4) = ssim_sub(NE_contact_mat,contact_mat,sub_window_size);
        error1(5) = ssim_sub(ND_contact_mat,contact_mat,sub_window_size);
        
        
        error2(1) = corr_sub(downsamped_contact_mat,contact_mat,sub_window_size);
        error2(2) = corr_sub(NR_F_contact_mat,contact_mat,sub_window_size);
        error2(3) = corr_sub(NR_F_contact_mat_twice,contact_mat,sub_window_size);
        error2(4) = corr_sub(NE_contact_mat,contact_mat,sub_window_size);
        error2(5) = corr_sub(ND_contact_mat,contact_mat,sub_window_size);
        
        
        error3(1) = corr_spearman_sub(downsamped_contact_mat,contact_mat,sub_window_size);
        error3(2) = corr_spearman_sub(NR_F_contact_mat,contact_mat,sub_window_size);
        error3(3) = corr_spearman_sub(NR_F_contact_mat_twice,contact_mat,sub_window_size);
        error3(4) = corr_spearman_sub(NE_contact_mat,contact_mat,sub_window_size);
        error3(5) = corr_spearman_sub(ND_contact_mat,contact_mat,sub_window_size);
        
        
    end
    
    error1 = error1_old + error1/loops;
    error2 = error2_old + error2/loops;
    error3 = error3_old + error3/loops;
end

toc




%% plot the performance
set(gcf,'color','white')
figure(1)
h1 = bar(error1);
ylabel('SSIM','FontSize',12)
title(['SSIM in GM12878 :',resolution],'FontSize',12)
set(gca,'XTickLabel',{'downsampled';'NR-F';'NR-F*2';'NE';'ND'})
set(gca,'XTickLabelRotation',90)
set(gca,'FontSize',12,'FontName','bold')
grid on;

set(gcf,'color','white')
figure(2)
h2 = bar(error2);
ylabel('Pearson correlation','FontSize',12)
title(['Pearson correlation in GM12878 :',resolution],'FontSize',12)
set(gca,'XTickLabel',{'downsampled';'NR-F';'NR-F*2';'NE';'ND'})
set(gca,'XTickLabelRotation',90)
set(gca,'FontSize',12,'FontName','bold')
grid on;

set(gcf,'color','white')
figure(3)
h3 = bar(error3);
ylabel('Spearman correlation','FontSize',12)
title(['Spearman correlation in GM12878 :',resolution],'FontSize',12)
set(gca,'XTickLabel',{'downsampled';'NR-F';'NR-F*2';'NE';'ND'})
set(gca,'XTickLabelRotation',90)
set(gca,'FontSize',12,'FontName','bold')
grid on;

