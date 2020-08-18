function c = MMSE_QR( H,x,snr )
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
    %��HH_e����qr�ֽ�
    [Q,R]=qr(HH);
    Q1=Q(1:NR,1:NT);
    Q2=Q(1:NT,1:NT);
    
    %��������������ź�
    y=Q1'*x(:,j);
    c(NT,j)=y(NT)/R(NT,NT);
    c(NT,j)=(c(NT,j)>=0)-(c(NT,j)<0)+0; 
    for k=NT-1:-1:1
        d=0;
        for i=k+1:NT
            d=d+R(k,i)*c(i,j);
        end
        z=y(k)-d;
        z=z/R(k,k);
        c(k,j)=(z>=0)-(z<0)+0; 
    end
end
c=(c+1)/2;
end
