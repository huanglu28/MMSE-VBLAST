function Q = swapc( Q,i,k ) 
%  交换Q的第i列和第k列
temp=Q(:,k);
Q(:,k)=Q(:,i);
Q(:,i)=temp;
end

