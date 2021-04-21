function score=corr_spearman_sub(contact_mat1,contact_mat2,sub_windows_size)

n=length(contact_mat1);
n_sub_region=n-sub_windows_size+1;
corr_vec=zeros(1,n_sub_region);

for i=1:n_sub_region
    region_left=i;
    region_right=region_left+sub_windows_size-1;
    
    mat1=contact_mat1(region_left:region_right,region_left:region_right);
    mat2=contact_mat2(region_left:region_right,region_left:region_right);
    
    corr_vec(i)=corr(mat1(:),mat2(:),'type' , 'Spearman');
end

score=mean(corr_vec(~isnan(corr_vec)));

end