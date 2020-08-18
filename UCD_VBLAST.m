%ģ��MIMOϵͳ 
%����������NT,����������NR��������󳤶�L
%NR>NT
%x=H*c+v
NT=4;
NR=4;
L=1000;
SNR=[0:3:30];%����ȣ�dB��
x_real=randint(NT,L);%NT*L�����ź�
x=zeros(NT,L);%��UCD�㷨�����õ��ķ����ź�

%ʵ�ʷ����źŵ�0ת��Ϊ-1,1����1
X=(-1).^(x_real+1);  

%%%%%%%%%%%%%%MIMO�ŵ�����
%�ŵ�ģ�ͣ�y=HFx+z
y=zeros(NR,L);

%���� H:��˥����NR*NT*Lά�����ŵ�
H=randn(NR,NT,L)+1i*randn(NR,NT,L);

%���Ӿ�ֵΪ0,����Ϊ1����̬�ֲ���NR*1ά�ĸ�˹������v
v=randn(NR,L)+1i*randn(NR,L);

erate=[];

for m=SNR
    snr=10^(m/10);
    z=1/snr*v;
    
    %%%%%����F����%%%%%
    for(l=1:L)
        
    %%%%%%%%%%%%%%%% transmitter��ģ�ⷢ�Ͷ� %%%%%%%%%%%%%%%
        %step1����H����svd�ֽ�
        [U,S,V]=svd(H(:,:,l));

        %step2�����ʷ������I(��λ�����)
        I=1;

        %step3��
        M=S*I;

        %step4:
        J=[U*M;sqrt(1/snr)*ones(NT,NT)];
        [~,M_J,~]=svd(J);
        M_J=M_J(1:NT,1:NT);
        
        %step5:��M_J���ξ�ֵ�ֽ�
        [Q,R,P]=gmdv(M_J);
        
        %step6:����F
        F=V*P;
        
        %���������ź�y
        y(:,l)=H(:,:,l)*F*X(:,l)+z(:,l);
    
    %%%%%%%%%%%%%%%% receiver:ģ����ն� %%%%%%%%%%%%%%%%
        %step7:����G=HF ��QR�ֽ�
        Q_GA=U*M*inv(M_J)*Q;
        %step8
        HH=Q_GA*R;
        G=inv(HH'*HH+(1/snr)*eye(NT))*HH';
        w=G*y(:,l);
        x(:,l)=(w>=0)-(w<0)+0;
    end
    x=(x+1)/2;
    %����UCD�㷨��������
    [errbit,err_ratio]=biterr(x_real,x);
    erate=[erate,err_ratio];
end
semilogy(SNR,erate,'d-b'); %��ɫ����
xlabel('SNR/dB');
ylabel('BER');
title('NT=4��NR=4ʱUCD�㷨�������ʺ�����ȹ�ϵ����');

