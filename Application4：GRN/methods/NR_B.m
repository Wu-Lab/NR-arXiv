function [output_network]=NETRW_B(mat,m)

[n_tf,n]=size(mat);
for i=1:n_tf
    mat(i,i)=0;
end

%% ********************** input matrix imputation ***************** 
mat(1:n_tf,1:n_tf)=(mat(1:n_tf,1:n_tf)+mat(1:n_tf,1:n_tf)')/2;
mat1=[mat;zeros(n-n_tf,n)]; 
mat1=(mat1+mat1')/2;
mat1=(mat1-min(mat1(:)))/(max(mat1(:))-min(mat1(:)));
mat1=mat1+min(mat1(mat1>0));

%% ********************** NETRW-B *********************************
P1=mat1./sum(mat1,2);
P2 = m * P1 /((m-1)*eye(n) + P1);
P2 = P2 - min(min(transpose(P2)),0)'; 
P2 = P2 ./ sum(P2,2);
stat_d=abs(null((P2-eye(n))'));
net_new=diag(stat_d)*P2;

%% ****************************************************************
net_new=net_new+net_new';
output_network=net_new/max(max(net_new));
output_network=output_network-diag(diag(output_network));
output_network=output_network(1:n_tf,:);

