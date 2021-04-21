function wn=NE_dn(w,type)
w = w*length(w);
w = double(w);
D = sum(abs(w),2)+eps;% 行和

if type == 'ave'
    D=1./D;
    D=sparse(1:length(D),1:length(D),D); %sparse和full是反义词
    wn=D*w;
elseif type == 'gph'
    D=1./sqrt(D);
    D=sparse(1:length(D),1:length(D),D);
    wn=D*(w*D);
end