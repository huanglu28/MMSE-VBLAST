function [ w,k ] = minnormc( Q,l,r )
%返回Q在l~r列中的最小范数列所在的列数k和该列的范数w
w=norm(Q(:,l),1);
k=l;
%计算Q的所有列范数
for i=l+1:r
    norm_q=norm(Q(:,i),1);
    if norm_q<w
        w=norm_q;
        k=i;
    end 
end
end

