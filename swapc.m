function Q = swapc( Q,i,k ) 
%  ����Q�ĵ�i�к͵�k��
temp=Q(:,k);
Q(:,k)=Q(:,i);
Q(:,i)=temp;
end

