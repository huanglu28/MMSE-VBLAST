%模拟MIMO系统 
%发射天线数NT,接收天线数NR，发射矩阵长度L
%NR>NT
%x=H*c+v
NT=4;
NR=6;
L=1000;
SNR=[0:1:20];%信噪比（dB）
c_real=randint(NT,L);%NT*L发射信号
c=zeros(NT,L);%经VBLAST算法计算后得到的发射信号


%实际发射信号的0转化为-1,1保持1
X=(-1).^(c_real+1);  

%%%%%%%%%%%%%%MIMO信道传输
%快衰弱的NR*NT*L维瑞利信道
H=sqrt(1/2)*(randn(NR,NT,L)+1i*randn(NR,NT,L));
%服从均值为0,方差为1的正态分布的NR*1维的高斯白噪声v
v=sqrt(1/2)*(randn(NR,L)+1i*randn(NR,L));

%未叠加噪声的接收信号x
x=zeros(NR,L);
for i=1:L
    x(:,i)=sqrt(1/2)*H(:,:,i)*X(:,i);
end

%%%%%%%%%%%%%%%%% V-blast算法 %%%%%%%%%%%%%%%%%
disp('V-blast算法');
%不同信噪比下的误码率
erate=[];
%叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %经解码得到的信号
    c=V_blast(H,x_noised);
    %计算V-blast算法的误码率
    [errbit,err_ratio]=biterr(c_real,c);
    erate=[erate,err_ratio];
end
semilogy(SNR,erate,'d-b'); %蓝色菱形
hold on;

%%%%%%%%%%%%%%%%% USQR 算法 %%%%%%%%%%%%%%%%%
disp('Unsorted QR deposition');
%不同信噪比下的误码率
erate_usqr=[];
%叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %经解码得到的信号
    c=USQR(H,x_noised);
    %计算USQR算法的误码率
    [errbit,err_ratio]=biterr(c_real,c);
    erate_usqr=[erate_usqr,err_ratio];
end
semilogy(SNR,erate_usqr,'--k'); %黑色星号
hold on;    

%%%%%%%%%%%%%%%%% SQRD算法 %%%%%%%%%%%%%%%%%
disp('Sorted QR deposition');
%不同信噪比下的误码率
erate_sqrd=[];
%叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %经解码得到的信号
    c=SQRD(H,x_noised);
    %计算SQRD算法的误码率
    [errbit,err_ratio]=biterr(c_real,c);
    erate_sqrd=[erate_sqrd,err_ratio];
end
semilogy(SNR,erate_sqrd,'o-g'); %绿色圆圈
hold on; 

%%%%%%%%%%%%%%%%% GMD算法 %%%%%%%%%%%%%%%%%
disp('Geometric Mean Decomposition');
%不同信噪比下的误码率
erate_gmd=[];
%叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %经解码得到的信号
    C=GMD(H,x_noised);
%计算SQRD算法的误码率
    [errbit,err_ratio]=biterr(c_real,C);
    erate_gmd=[erate_gmd,err_ratio];
end
semilogy(SNR,erate_gmd,'x-r'); %红色交叉
hold on; 

%%%%%%%%%%%%%%%%% MMSE算法 %%%%%%%%%%%%%%%%%
disp('MMSE算法');
%不同信噪比下的误码率
erate=[];
%叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %经解码得到的信号
    c=MMSE(H,x_noised,snr);
    %计算V-blast算法的误码率
    [errbit,err_ratio]=biterr(c_real,c);
    erate=[erate,err_ratio];
end
semilogy(SNR,erate,'d-r'); %红色菱形
hold on;

%%%%%%%%%%%%%%%%% MMSE_QR算法 %%%%%%%%%%%%%%%%%
disp('MMSE_QR算法');
%不同信噪比下的误码率
erate=[];
%叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %经解码得到的信号
    c=MMSE_QR(H,x_noised,snr);
    %计算MMSE_QR算法的误码率
    [errbit,err_ratio]=biterr(c_real,c);
    erate=[erate,err_ratio];
end
semilogy(SNR,erate,'--m'); %紫红色实形
hold on;

%%%%%%%%%%%%%%%%% MMSE_SQRD算法 %%%%%%%%%%%%%%%%%
disp('MMSE_SQRD算法');
%不同信噪比下的误码率
erate=[];
%叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %经解码得到的信号
    c=MMSE_SQRD(H,x_noised,snr);
    %计算MMSE_QR算法的误码率
    [errbit,err_ratio]=biterr(c_real,c);
    erate=[erate,err_ratio];
end
semilogy(SNR,erate,'*-k'); %紫红色实形
hold on;

xlabel('SNR');
ylabel('BER');
title('NT=4，NR=6时,MMSE算法与ZF算法的误码率和信噪比关系曲线');
legend('V-blast','unsorted QR','sorted QR','GMD','MMSE-BLAST','MMSE-QR','MMSE-SQRD');


   
           
        
        

