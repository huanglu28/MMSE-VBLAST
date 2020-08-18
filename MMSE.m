function c = MMSE( H,x,snr )
% ����MMSE׼��������㷨�����ʱ���ź�
% H -- NR*NTά�����ŵ�
% x -- �����ź�
% c -- �����ź�
% snr -- ��˹����������
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

