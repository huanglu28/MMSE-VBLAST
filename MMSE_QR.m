function c = MMSE_QR( H,x,snr )
% 无排序的MMSE算法
% H -- NR*NT维瑞利信道
% x -- 接收信号
% c -- 解码信号
% snr -- 高斯白噪声方差
[NR,NT,L]=size(H);
c=zeros(NT,L);
for j=1:L
    HH=H(:,:,j);
    %HH_e -- extended channel matrix
    HH=[HH;sqrt(1/snr)*ones(NT,NT)];
    %对HH_e进行qr分解
    [Q,R]=qr(HH);
    Q1=Q(1:NR,1:NT);
    Q2=Q(1:NT,1:NT);
    
    %后向带入解出发送信号
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
