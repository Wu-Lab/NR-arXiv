function [ W_out] = Network_Enhancement( W_in,order,K,alpha)
%%% W_in the the input network of size N X N
%%% K the number of neighbors
%%% alpha is the regularization parameter
%%% order determines the extent of diffusion. Typical values are 0.5,1,2.


if nargin<2
    K = min(20,ceil(length(W_in)/10));
    alpha = 0.9;
    order = 2;
end

if nargin<3
    K =min(20,ceil(length(W_in)/10));
    alpha = 0.9;
end

if nargin<4
    alpha = 0.9;
end


zeroindex = find(sum(abs(W_in))>0);
W0 = W_in(zeroindex,zeroindex);% 取正
W = NE_dn(W0,'ave');% 行归一化
W = (W+W')/2;

DD = sum(abs(W0));

if length(unique(W(:)))==2% 矩阵W中只有两个值？
    P = W;
else
    P = (dominateset(double(abs(W)),min(K,length(W)-1))).*sign(W);
    %只保留前K个邻居值并对称化
end


P = P + (eye(length(P))+diag(sum(abs(P'))));
%P = P + (eye(length(P)));
P = TransitionFields(P);

%P即为W的局部矩阵，即论文中的T矩阵
[U,D] = eig(full(P));
d = real(diag(D-eps));% 特征值
d = (1-alpha)*d./(1-alpha*d.^order);% 特征值变换，order是阶数
D = diag(real(d));
W = U*D*U';

W = W./repmat(1-diag(W),1,length(W));
D=sparse(1:length(DD),1:length(DD),(DD));
W=D*(W);
W(W<0)=0;

W = (W+W')/2;
W_out = zeros(size(W_in));
W_out(zeroindex,zeroindex) = W;
end

