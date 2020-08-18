function c = MMSE( H,x,snr )
% 基于MMSE准则的线性算法检测多层时空信号
% H -- NR*NT维瑞利信道
% x -- 接收信号
% c -- 解码信号
% snr -- 高斯白噪声方差
[NR,NT,L]=size(H);
c=zeros(NT,L);
for j=1:L
    HH=H(:,:,j);
    %G -- filter matrix
    G=inv(HH'*HH+(1/snr)*eye(NT))*HH';
    y=G*x(:,j);
    c(:,j)=(y>=0)-(y<0)+0;
end
c=(c+1)/2;
end

