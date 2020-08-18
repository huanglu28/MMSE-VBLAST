%模拟MIMO系统 
%发射天线数NT,接收天线数NR，发射矩阵长度L
%NR>NT
%x=H*c+v
NT=4;
NR=4;
L=1000;
SNR=[0:3:30];%信噪比（dB）
x_real=randint(NT,L);%NT*L发射信号
x=zeros(NT,L);%经UCD算法计算后得到的发射信号

%实际发射信号的0转化为-1,1保持1
X=(-1).^(x_real+1);  

%%%%%%%%%%%%%%MIMO信道传输
%信道模型：y=HFx+z
y=zeros(NR,L);

%生成 H:快衰弱的NR*NT*L维瑞利信道
H=randn(NR,NT,L)+1i*randn(NR,NT,L);

%服从均值为0,方差为1的正态分布的NR*1维的高斯白噪声v
v=randn(NR,L)+1i*randn(NR,L);

erate=[];

for m=SNR
    snr=10^(m/10);
    z=1/snr*v;
    
    %%%%%生成F矩阵%%%%%
    for(l=1:L)
        
    %%%%%%%%%%%%%%%% transmitter：模拟发送端 %%%%%%%%%%%%%%%
        %step1：对H进行svd分解
        [U,S,V]=svd(H(:,:,l));

        %step2：功率分配矩阵I(单位阵代替)
        I=1;

        %step3：
        M=S*I;

        %step4:
        J=[U*M;sqrt(1/snr)*ones(NT,NT)];
        [~,M_J,~]=svd(J);
        M_J=M_J(1:NT,1:NT);
        
        %step5:对M_J几何均值分解
        [Q,R,P]=gmdv(M_J);
        
        %step6:产生F
        F=V*P;
        
        %产生接收信号y
        y(:,l)=H(:,:,l)*F*X(:,l)+z(:,l);
    
    %%%%%%%%%%%%%%%% receiver:模拟接收端 %%%%%%%%%%%%%%%%
        %step7:计算G=HF 的QR分解
        Q_GA=U*M*inv(M_J)*Q;
        %step8
        HH=Q_GA*R;
        G=inv(HH'*HH+(1/snr)*eye(NT))*HH';
        w=G*y(:,l);
        x(:,l)=(w>=0)-(w<0)+0;
    end
    x=(x+1)/2;
    %计算UCD算法的误码率
    [errbit,err_ratio]=biterr(x_real,x);
    erate=[erate,err_ratio];
end
semilogy(SNR,erate,'d-b'); %蓝色菱形
xlabel('SNR/dB');
ylabel('BER');
title('NT=4，NR=4时UCD算法的误码率和信噪比关系曲线');

