function [W_output] = NR_F(W_input,m)

n=length(W_input);

P1=W_input./sum(W_input,2);
P1=fillmissing(P1,'constant',0);
P2=(m-1)*P1/(m*eye(n)-P1);

stat=null((P2-eye(n))');
check=size(stat);

% if the stationary distribution dont unique
% we use the Random Walk with Restart (RWR) instead
if check(2)~=1
    P2=P2+1e-10;
    P2=P2./sum(P2,2);
    stat=null((P2-eye(n))');
end

% double check
check=size(stat);
if check(2)~=1
    for i=1:check(2)
        if sum(stat(:,i)>=0)==n||sum(stat(:,i)<=0)==n
            keep_index=i;
        end
    end
    stat=stat(:,keep_index);
end

W_output=diag(abs(stat))*P2;
W_output=W_output*max(max(W_input))/max(max(W_output));

end

