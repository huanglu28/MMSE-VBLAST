function [ w,k ] = minnormc( Q,l,r )
%����Q��l~r���е���С���������ڵ�����k�͸��еķ���w
w=norm(Q(:,l),1);
k=l;
%����Q�������з���
for i=l+1:r
    norm_q=norm(Q(:,i),1);
    if norm_q<w
        w=norm_q;
        k=i;
    end 
end
end

