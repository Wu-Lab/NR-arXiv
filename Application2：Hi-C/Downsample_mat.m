function contact_mat_downsamped=Downsample_mat(reads_mat,downsample_ratio)

[n1,n2]=size(reads_mat);
contact_mat_downsamped=zeros(n1,n2);

total_reads=sum(sum(reads_mat));
num_downsampled = round(total_reads*downsample_ratio);

sample_index=randsample(total_reads,num_downsampled);

ifsampled=zeros(1,total_reads);
ifsampled(sample_index)=1;


bin_left=1;

for i=1:n1
    for j=1:n2
        if reads_mat(i,j)~=0
            bin_right=reads_mat(i,j) + bin_left;
            contact_mat_downsamped(i,j)=sum(sum(ifsampled(bin_left:(bin_right-1))));
            bin_left=bin_right;
        end
    end
end
        
end