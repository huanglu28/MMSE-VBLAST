function C = GMD( H,x )
% ���ξ�ֵ�ֽ����MIMOϵͳ�еķ����ź�
% H -- NR*NTά�����ŵ�
% x -- �����ź�
% c -- �����ź�
[NR,NT,L]=size(H);
c=zeros(NT,L);
    C=zeros(NT,L);
for j=1:L
    HH=H(:,:,j);
%     
%     %��H����SVD�ֽ�
%     [U,S,V]=svd(HH);
%     %R�ĶԽ���Ԫ�ص���HH������ֵ
%     R=S;
%     %ȡ��R�ĶԽ���Ԫ��
%     d=diag(R);
%     %����H������ֵrH
%     rH=1;
%     for i=1:K
%         rH=rH*abs(d(i));
%     end
%     rH=power(rH,1/K);
%     %����Q��P
%     Q=U;
%     P=V;
%     z=zeros(K-1,1);
%     for m=1:K-1
%         %R�ĵ�m���Խ���Ԫ�ؼ�Ϊa1��m+1����Ϊa2
%         a1=d(m);
%         a2=d(m+1);
%         if a1==a2&&a1==rH
%             c=1;
%             s=0;
%         else
%             c=sqrt((rH^2-a1^2)/(a1^2-a2^2));
%             s=sqrt(1-c^2);
%         end
%         d(m+1)=a1*a2/rH;  % y
%         z(m)=s*c*(a2^2-a1^2)/rH;
%         %��R�ĶԽ���Ԫ�ػ�ΪrH
%         R(m,m)=rH;
%         R(1:m-1,m)=z(1:m-1)*c;
%         z(1:m-1)=-z(1:m-1)*s;
%          %���ɼ���˹��ת����G1��G2
%         G1=[c,-s;s,c];
%         G2=1/rH*[c*a1,s*a2;-s*a2,c*a1];
%         Q(:,[m m+1])=Q(:,[m m+1])*G2';
%         P(:,[m m+1])=P(:,[m m+1])*G1;
%     end
%     R(K,K)=rH;
%     R(1:K-1,K)=z;
%     %x=Hc+z
%     %x=QRP'*c+z
%     %y=Q'x,c=Ps
%     %y=Rs+Q'z
%     %��s
    K=NT;
    [Q,R,P]=gmdv(HH);
    y=Q'*x(:,j);
    S=zeros(K,1);
    S(K)=y(K)/(sqrt(1/2)*R(K,K));
    for n=K-1:-1:1
        sum=0;
        for o=n+1:K
            sum=sum+R(n,o)*S(o)*sqrt(1/2);
        end
        S(n)=(y(n)-sum)/R(n,n);
    end
    C(:,j)=P*S;
end
C=(C>=0)-(C<0)+0;
C=(C+1)/2;
end
