%模拟MIMO系统 
%发射天线数NT,接收天线数NR，发射矩阵长度L
%NR>NT
%x=H*c+v
NT=4;
NR=4;
L=1000;
SNR=[0:1:30];%信噪比（dB）
c_real=randint(NT,L);%NT*L发射信号
c=zeros(NT,L);%经VBLAST算法计算后得到的发射信号


%实际发射信号的0转化为-1,1保持1
X=(-1).^(c_real+1);  

%%%%%%%%%%%%%%MIMO信道传输
%快衰弱的NR*NT*L维瑞利信道
H=randn(NR,NT,L)+1i*randn(NR,NT,L);
%服从均值为0,方差为1的正态分布的NR*1维的高斯白噪声v
v=randn(NR,L)+1i*randn(NR,L);

%未叠加噪声的接收信号x
x=zeros(NR,L);
for i=1:L
    x(:,i)=H(:,:,i)*X(:,i);
end
% 
% %%%%%%%%%%%%%%%%% MMSE算法 %%%%%%%%%%%%%%%%%
% disp('MMSE算法');
% %不同信噪比下的误码率
% erate=[];
% %叠加噪声
% for m=SNR
%     snr=10^(m/10);
%     x_noised=x+sqrt(1/snr)*v;
%     %经解码得到的信号
%     c=MMSE(H,x_noised,snr);
%     %计算V-blast算法的误码率
%     [errbit,err_ratio]=biterr(c_real,c);
%     erate=[erate,err_ratio];
% end
% semilogy(SNR,erate,'d-r'); %红色菱形
% hold on;
% 
% %%%%%%%%%%%%%%%%% MMSE_QR算法 %%%%%%%%%%%%%%%%%
% disp('MMSE_QR算法');
% %不同信噪比下的误码率
% erate=[];
% %叠加噪声
% for m=SNR
%     snr=10^(m/10);
%     x_noised=x+sqrt(1/snr)*v;
%     %经解码得到的信号
%     c=MMSE_QR(H,x_noised,snr);
%     %计算MMSE_QR算法的误码率
%     [errbit,err_ratio]=biterr(c_real,c);
%     erate=[erate,err_ratio];
% end
% semilogy(SNR,erate,'--m'); %紫红色实形
% hold on;
% 
% %%%%%%%%%%%%%%%%% MMSE_SQRD算法 %%%%%%%%%%%%%%%%%
% disp('MMSE_SQRD算法');
% %不同信噪比下的误码率
% erate=[];
% %叠加噪声
% for m=SNR
%     snr=10^(m/10);
%     x_noised=x+sqrt(1/snr)*v;
%     %经解码得到的信号
%     c=MMSE_SQRD(H,x_noised,snr);
%     %计算MMSE_QR算法的误码率
%     [errbit,err_ratio]=biterr(c_real,c);
%     erate=[erate,err_ratio];
% end
% semilogy(SNR,erate,'*-k'); %紫红色实形
% hold on;

%%%%%%%%%%%%%%%%% MMSE_SQRD_PSA算法 %%%%%%%%%%%%%%%%%
disp('MMSE_SQRD_PSA算法');
%不同信噪比下的误码率
erate=[];
%叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %经解码得到的信号
    c=MMSE_SQRD_PSA(H,x_noised,snr);
    %计算MMSE_QR算法的误码率
    [errbit,err_ratio]=biterr(c_real,c);
    erate=[erate,err_ratio];
end
semilogy(SNR,erate,'k'); %黑色直线
hold on;

xlabel('SNR/dB');
ylabel('BER');
title('NT=4，NR=4时三种MMSE算法的误码率和信噪比关系曲线');
legend('MMSE','MMSE-QR','MMSE-SQRD');

