function C = GMD( H,x )
% 几何均值分解解码MIMO系统中的发射信号
% H -- NR*NT维瑞利信道
% x -- 接收信号
% c -- 解码信号
[NR,NT,L]=size(H);
c=zeros(NT,L);
    C=zeros(NT,L);
for j=1:L
    HH=H(:,:,j);
%     
%     %对H进行SVD分解
%     [U,S,V]=svd(HH);
%     %R的对角线元素等于HH的特征值
%     R=S;
%     %取出R的对角线元素
%     d=diag(R);
%     %计算H的奇异值rH
%     rH=1;
%     for i=1:K
%         rH=rH*abs(d(i));
%     end
%     rH=power(rH,1/K);
%     %生成Q、P
%     Q=U;
%     P=V;
%     z=zeros(K-1,1);
%     for m=1:K-1
%         %R的第m个对角线元素记为a1，m+1个记为a2
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
%         %将R的对角线元素化为rH
%         R(m,m)=rH;
%         R(1:m-1,m)=z(1:m-1)*c;
%         z(1:m-1)=-z(1:m-1)*s;
%          %生成吉文斯旋转矩阵G1、G2
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
%     %求s
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
