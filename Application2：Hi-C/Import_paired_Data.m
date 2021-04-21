function [contact_mat,downsamped_contact_mat] = Import_paired_Data(chromosome,resolution,downsample_ratio)

disp(['******loading data with ',resolution,' resolution in chromosome ',chromosome(isstrprop(chromosome,'digit')),'*****']);

%% load path
reads_path=['./example data/',chromosome,'_',resolution,'.RAWobserved'];
nor_path=['./example data/',chromosome,'_',resolution,'.KRnorm'];

% reads_path=[path_par,'/',resolution,'_resolution_intrachromosomal/', ...
%     chromosome,'/MAPQGE30/',chromosome,'_',resolution,'.RAWobserved'];
% nor_path=[path_par,'/',resolution,'_resolution_intrachromosomal/', ...
%     chromosome,'/MAPQGE30/',chromosome,'_',resolution,'.KRnorm'];

%% load data

reads_data_3col=load(reads_path);
nor_data=load(nor_path);
nor_data(nor_data==0)=1e-10;
n_bin=length(nor_data);

nan_index=find(isnan(nor_data));

%% convert into raw data matrix
%disp('convert into matrix.....')
resolution_num=str2num(resolution(isstrprop(resolution,'digit')))*1000;
n_contact=length(reads_data_3col(:,1));  %number of contacts between different locus
reads_mat=zeros(n_bin,n_bin);

for i=1:n_contact
    index_row=reads_data_3col(i,1)/resolution_num+1;
    index_col=reads_data_3col(i,2)/resolution_num+1;
    contact=reads_data_3col(i,3);
    reads_mat(index_row,index_col)=contact;
end

% ÊÍ·ÅÄÚ´æ
clear reads_data_3col

reads_mat=reads_mat+reads_mat'-diag(diag(reads_mat));


%% downsample
downsamped_reads_mat=sub_Downsample_mat(reads_mat,downsample_ratio);
downsamped_reads_mat=(downsamped_reads_mat+downsamped_reads_mat')/2;


%% keep valid value
nor_data(nan_index)=[];
reads_mat(nan_index,:)=[];
reads_mat(:,nan_index)=[];
downsamped_reads_mat(nan_index,:)=[];
downsamped_reads_mat(:,nan_index)=[];


%% normalize

downsamped_contact_mat=downsamped_reads_mat./(nor_data*nor_data');
contact_mat=reads_mat./(nor_data*nor_data');

%%
downsamped_contact_mat(downsamped_contact_mat<1e-6)=0;
contact_mat(contact_mat<1e-6)=0;

zero_index=find(sum(downsamped_contact_mat)==0);
downsamped_contact_mat(zero_index,:)=[];
downsamped_contact_mat(:,zero_index)=[];
contact_mat(zero_index,:)=[];
contact_mat(:,zero_index)=[];

disp('********************Loading completed************************')
end