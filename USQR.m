function c = USQR( H,x )
% unsorted QR�㷨�����ʱ���ź�
% H -- NR*NTά�����ŵ�
% x -- �����ź�
% c -- �����ź�
[NR,NT,L]=size(H);
c=zeros(NT,L);
for j=1:L
    R=zeros(NT,NT);
    Q=H(:,:,j);
    for i=1:NT
        %R�ĶԽ���Ԫ�ص���qi��1����
        R(i,i)=norm(Q(:,i),1);
        Q(:,i)=Q(:,i)/R(i,i);
        for l=i+1:NT
            R(i,l)=Q(:,i)'*Q(:,l);
            Q(:,l)=Q(:,l)-R(i,l)*Q(:,i);
        end
    end
    
    y=Q'*x(:,j);
    c(NT,j)=y(NT)/R(NT,NT);
    c(NT,j)=(c(NT,j)>=0)-(c(NT,j)<0)+0; 
    for k=NT-1:-1:1
        d=0;
        for i=k+1:NT
            d=d+sqrt(1/NT)*R(k,i)*c(i,j);
        end
        z=y(k)-d;
        z=z/R(k,k);
        c(k,j)=(z>=0)-(z<0)+0;  
    end
end
c=(c+1)/2;
end



