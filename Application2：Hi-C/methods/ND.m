function mat_nd=ND(contact_mat)

n = size(contact_mat,1);
contact_mat = (contact_mat+contact_mat')/2;

%***********************************
% eigen decomposition
[U,D] = eig(contact_mat);

lam_n=abs(min(min(diag(D)),0));
lam_p=abs(max(max(diag(D)),0));

beta = 0.5;
m1=lam_p*(1-beta)/beta;
m2=lam_n*(1+beta)/beta;
m=max(m1,m2);

%network deconvolution
for i = 1:n
    D(i,i) = (D(i,i))/(m+D(i,i));
end
mat_nd = U*D*inv(U);
m1 = min(min(mat_nd));
m2 = max(max(mat_nd));
mat_nd = (mat_nd-m1)/(m2-m1);






