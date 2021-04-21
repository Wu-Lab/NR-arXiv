function denoised_contact_mat_out=denoise(downsamped_contact_mat,sub_window_size,m,method_name)

gap=max(2,round(sub_window_size/20));
n=length(downsamped_contact_mat);
from_max=n-sub_window_size+1; %区间左侧能取到的最大值
region_size=length(1:gap:from_max);
denoised_contact_mat=zeros(n,n,region_size);

for i=1:region_size
    from=(i-1)*gap+1;
    to=from+sub_window_size-1;
        
        if method_name==1
            denoised_contact_mat(from:to,from:to,i)=...
                NR_F(downsamped_contact_mat(from:to,from:to),m);
        elseif method_name==2
            denoised_contact_mat(from:to,from:to,i)=...
                Network_Enhancement(downsamped_contact_mat(from:to,from:to));
        elseif method_name==3
            denoised_contact_mat(from:to,from:to,i)=...
                ND(downsamped_contact_mat(from:to,from:to));
        else
            disp('Wrong method!');
        end
    
end

denoised_contact_mat_out=zeros(n,n);
for i=1:n
    for j=1:n
        temp=denoised_contact_mat(i,j,:);
        denoised_contact_mat_out(i,j)=mean(temp(temp>0));
    end
end

denoised_contact_mat_out=fillmissing(denoised_contact_mat_out,'constant',0);

end