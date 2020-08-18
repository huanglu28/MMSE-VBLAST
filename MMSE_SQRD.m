function c = MMSE_SQRD( H,x,snr )
% �������MMSE�㷨
% H -- NR*NTά�����ŵ�
% x -- �����ź�
% c -- �����ź�
% snr -- ��˹����������
[NR,NT,L]=size(H);
c=zeros(NT,L);
for j=1:L
    HH=H(:,:,j);
    %HH_e -- extended channel matrix
    HH=[HH;sqrt(1/snr)*ones(NT,NT)];
    
    %��HH����������qr�ֽ�
    R=zeros(NT,NT);
    Q=HH;
    S=1:NT;
    for i=1:NT
    %�ҳ�Q����С�����к͸��еķ���
    [min_norm,k]=minnormc(Q,i,NT);
    
        if(i~=k)
        %����Q��R��S�ĵ�i�к͵�k��
            Q=swapc(Q,i,k);
            R=swapc(R,i,k);
            S=swapc(S,i,k);
        end
       
        %R�ĶԽ���Ԫ�ص�����С������1����
        R(i,i)=min_norm;
        Q(:,i)=Q(:,i)/R(i,i);
        for l=i+1:NT
            R(i,l)=Q(:,i)'*Q(:,l);
            Q(:,l)=Q(:,l)-R(i,l)*Q(:,i);
        end
    end
    Q1=Q(1:NR,1:NT);
    Q2=Q(1:NT,1:NT);
   
    y=Q1'*x(:,j);
    c_u=zeros(NT,1);  %δ�û���c����
    %�����NT���źŵĴ�С
    c_u(NT)=y(NT)/R(NT,NT);
    c_u(NT)=(c_u(NT)>=0)-(c_u(NT)<0)+0; 

    for k=NT-1:-1:1
        d=0;
        for i=k+1:NT
            d=d+R(k,i)*c_u(i);
        end
        z=y(k)-d;
        z=z/R(k,k);
        c_u(k)=(z>=0)-(z<0)+0;  
    end
    
    for h=1:NT
        c(S(h),j)=c_u(h);
    end  
end
c=(c+1)/2;
end


