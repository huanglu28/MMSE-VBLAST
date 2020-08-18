%ģ��MIMOϵͳ 
%����������NT,����������NR��������󳤶�L
%NR>NT
%x=H*c+v
NT=4;
NR=6;
L=1000;
SNR=[0:1:20];%����ȣ�dB��
c_real=randint(NT,L);%NT*L�����ź�
c=zeros(NT,L);%��VBLAST�㷨�����õ��ķ����ź�


%ʵ�ʷ����źŵ�0ת��Ϊ-1,1����1
X=(-1).^(c_real+1);  

%%%%%%%%%%%%%%MIMO�ŵ�����
%��˥����NR*NT*Lά�����ŵ�
H=sqrt(1/2)*(randn(NR,NT,L)+1i*randn(NR,NT,L));
%���Ӿ�ֵΪ0,����Ϊ1����̬�ֲ���NR*1ά�ĸ�˹������v
v=sqrt(1/2)*(randn(NR,L)+1i*randn(NR,L));

%δ���������Ľ����ź�x
x=zeros(NR,L);
for i=1:L
    x(:,i)=sqrt(1/2)*H(:,:,i)*X(:,i);
end

%%%%%%%%%%%%%%%%% V-blast�㷨 %%%%%%%%%%%%%%%%%
disp('V-blast�㷨');
%��ͬ������µ�������
erate=[];
%��������
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %������õ����ź�
    c=V_blast(H,x_noised);
    %����V-blast�㷨��������
    [errbit,err_ratio]=biterr(c_real,c);
    erate=[erate,err_ratio];
end
semilogy(SNR,erate,'d-b'); %��ɫ����
hold on;

%%%%%%%%%%%%%%%%% USQR �㷨 %%%%%%%%%%%%%%%%%
disp('Unsorted QR deposition');
%��ͬ������µ�������
erate_usqr=[];
%��������
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %������õ����ź�
    c=USQR(H,x_noised);
    %����USQR�㷨��������
    [errbit,err_ratio]=biterr(c_real,c);
    erate_usqr=[erate_usqr,err_ratio];
end
semilogy(SNR,erate_usqr,'--k'); %��ɫ�Ǻ�
hold on;    

%%%%%%%%%%%%%%%%% SQRD�㷨 %%%%%%%%%%%%%%%%%
disp('Sorted QR deposition');
%��ͬ������µ�������
erate_sqrd=[];
%��������
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %������õ����ź�
    c=SQRD(H,x_noised);
    %����SQRD�㷨��������
    [errbit,err_ratio]=biterr(c_real,c);
    erate_sqrd=[erate_sqrd,err_ratio];
end
semilogy(SNR,erate_sqrd,'o-g'); %��ɫԲȦ
hold on; 

%%%%%%%%%%%%%%%%% GMD�㷨 %%%%%%%%%%%%%%%%%
disp('Geometric Mean Decomposition');
%��ͬ������µ�������
erate_gmd=[];
%��������
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %������õ����ź�
    C=GMD(H,x_noised);
%����SQRD�㷨��������
    [errbit,err_ratio]=biterr(c_real,C);
    erate_gmd=[erate_gmd,err_ratio];
end
semilogy(SNR,erate_gmd,'x-r'); %��ɫ����
hold on; 

%%%%%%%%%%%%%%%%% MMSE�㷨 %%%%%%%%%%%%%%%%%
disp('MMSE�㷨');
%��ͬ������µ�������
erate=[];
%��������
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %������õ����ź�
    c=MMSE(H,x_noised,snr);
    %����V-blast�㷨��������
    [errbit,err_ratio]=biterr(c_real,c);
    erate=[erate,err_ratio];
end
semilogy(SNR,erate,'d-r'); %��ɫ����
hold on;

%%%%%%%%%%%%%%%%% MMSE_QR�㷨 %%%%%%%%%%%%%%%%%
disp('MMSE_QR�㷨');
%��ͬ������µ�������
erate=[];
%��������
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %������õ����ź�
    c=MMSE_QR(H,x_noised,snr);
    %����MMSE_QR�㷨��������
    [errbit,err_ratio]=biterr(c_real,c);
    erate=[erate,err_ratio];
end
semilogy(SNR,erate,'--m'); %�Ϻ�ɫʵ��
hold on;

%%%%%%%%%%%%%%%%% MMSE_SQRD�㷨 %%%%%%%%%%%%%%%%%
disp('MMSE_SQRD�㷨');
%��ͬ������µ�������
erate=[];
%��������
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %������õ����ź�
    c=MMSE_SQRD(H,x_noised,snr);
    %����MMSE_QR�㷨��������
    [errbit,err_ratio]=biterr(c_real,c);
    erate=[erate,err_ratio];
end
semilogy(SNR,erate,'*-k'); %�Ϻ�ɫʵ��
hold on;

xlabel('SNR');
ylabel('BER');
title('NT=4��NR=6ʱ,MMSE�㷨��ZF�㷨�������ʺ�����ȹ�ϵ����');
legend('V-blast','unsorted QR','sorted QR','GMD','MMSE-BLAST','MMSE-QR','MMSE-SQRD');


   
           
        
        

