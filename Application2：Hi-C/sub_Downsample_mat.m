function contact_mat_downsamped=sub_Downsample_mat(reads_mat,downsample_ratio)

n=length(reads_mat);
contact_mat_downsamped=zeros(n,n);
sub_window_size=1000;

if n<sub_window_size
    contact_mat_downsamped=Downsample_mat(reads_mat,downsample_ratio);
else
    section=floor(n/sub_window_size);
    for i=1:section
        row_start=(i-1)*sub_window_size+1;
        row_end=row_start+sub_window_size-1;
        for j=1:section
            col_start=(j-1)*sub_window_size+1;
            col_end=col_start+sub_window_size-1;
            contact_mat_downsamped(row_start:row_end,col_start:col_end)=...
                Downsample_mat(reads_mat(row_start:row_end,col_start:col_end),downsample_ratio);
        end
    end
    
    if row_end~=n
        contact_mat_downsamped((row_end+1):n,1:n)=...
            Downsample_mat(reads_mat((row_end+1):n,1:n),downsample_ratio);
        contact_mat_downsamped(1:row_end,(col_end+1):n)=...
            Downsample_mat(reads_mat(1:row_end,(col_end+1):n),downsample_ratio);
        
    end
end


end