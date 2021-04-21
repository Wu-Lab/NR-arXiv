function [ Wnew ] = TransitionFields( W )
zeroindex = find(sum(W,2)==0);%行和等于0的
W = W*length(W);
W = NE_dn(W,'ave');%行归一化

w = sqrt(sum(abs(W))+eps);%列和开方
W = W./repmat(w,length(W),1);%除以列和

W = W*W';
Wnew = W;
Wnew(zeroindex,:) = 0;
Wnew(:,zeroindex) = 0;
end

