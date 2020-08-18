%ģ��MIMOϵͳ 
%����������NT,����������NR��������󳤶�L
%NR>NT
%x=H*c+v
NT=4;
NR=4;
L=1000;
SNR=[0:1:30];%����ȣ�dB��
c_real=randint(NT,L);%NT*L�����ź�
c=zeros(NT,L);%��VBLAST�㷨�����õ��ķ����ź�


%ʵ�ʷ����źŵ�0ת��Ϊ-1,1����1
X=(-1).^(c_real+1);  

%%%%%%%%%%%%%%MIMO�ŵ�����
%��˥����NR*NT*Lά�����ŵ�
H=randn(NR,NT,L)+1i*randn(NR,NT,L);
%���Ӿ�ֵΪ0,����Ϊ1����̬�ֲ���NR*1ά�ĸ�˹������v
v=randn(NR,L)+1i*randn(NR,L);

%δ���������Ľ����ź�x
x=zeros(NR,L);
for i=1:L
    x(:,i)=H(:,:,i)*X(:,i);
end
% 
% %%%%%%%%%%%%%%%%% MMSE�㷨 %%%%%%%%%%%%%%%%%
% disp('MMSE�㷨');
% %��ͬ������µ�������
% erate=[];
% %��������
% for m=SNR
%     snr=10^(m/10);
%     x_noised=x+sqrt(1/snr)*v;
%     %������õ����ź�
%     c=MMSE(H,x_noised,snr);
%     %����V-blast�㷨��������
%     [errbit,err_ratio]=biterr(c_real,c);
%     erate=[erate,err_ratio];
% end
% semilogy(SNR,erate,'d-r'); %��ɫ����
% hold on;
% 
% %%%%%%%%%%%%%%%%% MMSE_QR�㷨 %%%%%%%%%%%%%%%%%
% disp('MMSE_QR�㷨');
% %��ͬ������µ�������
% erate=[];
% %��������
% for m=SNR
%     snr=10^(m/10);
%     x_noised=x+sqrt(1/snr)*v;
%     %������õ����ź�
%     c=MMSE_QR(H,x_noised,snr);
%     %����MMSE_QR�㷨��������
%     [errbit,err_ratio]=biterr(c_real,c);
%     erate=[erate,err_ratio];
% end
% semilogy(SNR,erate,'--m'); %�Ϻ�ɫʵ��
% hold on;
% 
% %%%%%%%%%%%%%%%%% MMSE_SQRD�㷨 %%%%%%%%%%%%%%%%%
% disp('MMSE_SQRD�㷨');
% %��ͬ������µ�������
% erate=[];
% %��������
% for m=SNR
%     snr=10^(m/10);
%     x_noised=x+sqrt(1/snr)*v;
%     %������õ����ź�
%     c=MMSE_SQRD(H,x_noised,snr);
%     %����MMSE_QR�㷨��������
%     [errbit,err_ratio]=biterr(c_real,c);
%     erate=[erate,err_ratio];
% end
% semilogy(SNR,erate,'*-k'); %�Ϻ�ɫʵ��
% hold on;

%%%%%%%%%%%%%%%%% MMSE_SQRD_PSA�㷨 %%%%%%%%%%%%%%%%%
disp('MMSE_SQRD_PSA�㷨');
%��ͬ������µ�������
erate=[];
%��������
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %������õ����ź�
    c=MMSE_SQRD_PSA(H,x_noised,snr);
    %����MMSE_QR�㷨��������
    [errbit,err_ratio]=biterr(c_real,c);
    erate=[erate,err_ratio];
end
semilogy(SNR,erate,'k'); %��ɫֱ��
hold on;

xlabel('SNR/dB');
ylabel('BER');
title('NT=4��NR=4ʱ����MMSE�㷨�������ʺ�����ȹ�ϵ����');
legend('MMSE','MMSE-QR','MMSE-SQRD');

