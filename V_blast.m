function c =V_blast( H,x )
% V-blast�㷨�����ʱ���ź�
% H -- NR*NTά�����ŵ�
% x -- �����ź�
% c -- �����ź�
[NR,NT,L]=size(H);
c=zeros(NT,L);
    for j=1:L
        HH=H(:,:,j);
        S=1:NR;%��Ǳ���ȥ����
        for i=1:NT
            G=pinv(HH);%GΪH�Ķ�ά���������
            [w,k]=minnorm(G);%kΪG����С�з���������
            y=G(k,:)*x(:,j); %�����о�ͳ����
            temp=S(k);
            %��y�����о����õ������ź�
            c(temp,j)=1*(y>=0)-1*(y<0)+0;
            %�ӽ����ź�������c���õ���������ź�
            x(:,j)=x(:,j)-sqrt(1/NR)*HH(:,k)*c(temp,j);
            HH(:,k)=[];
            S(k)=[];
        end
        end
   c=(c+1)/2;
end






