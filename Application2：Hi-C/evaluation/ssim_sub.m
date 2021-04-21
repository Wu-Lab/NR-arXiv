function score=ssim_sub(contact_mat1,contact_mat2,sub_windows_size)

n=length(contact_mat1);
n_sub_region=n-sub_windows_size+1;
ssim_vec=zeros(1,n_sub_region);

for i=1:n_sub_region
    region_left=i;
    region_right=region_left+sub_windows_size-1;
    
    ssim_vec(i)=ssim(contact_mat1(region_left:region_right,region_left:region_right),...
        contact_mat2(region_left:region_right,region_left:region_right));
end


score=mean(ssim_vec);

end