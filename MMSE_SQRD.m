function c = MMSE_SQRD( H,x,snr )
% 有排序的MMSE算法
% H -- NR*NT维瑞利信道
% x -- 接收信号
% c -- 解码信号
% snr -- 高斯白噪声方差
[NR,NT,L]=size(H);
c=zeros(NT,L);
for j=1:L
    HH=H(:,:,j);
    %HH_e -- extended channel matrix
    HH=[HH;sqrt(1/snr)*ones(NT,NT)];
    
    %对HH进行有排序qr分解
    R=zeros(NT,NT);
    Q=HH;
    S=1:NT;
    for i=1:NT
    %找出Q的最小范数列和该列的范数
    [min_norm,k]=minnormc(Q,i,NT);
    
        if(i~=k)
        %交换Q、R、S的第i列和第k列
            Q=swapc(Q,i,k);
            R=swapc(R,i,k);
            S=swapc(S,i,k);
        end
       
        %R的对角线元素等于最小范数列1范数
        R(i,i)=min_norm;
        Q(:,i)=Q(:,i)/R(i,i);
        for l=i+1:NT
            R(i,l)=Q(:,i)'*Q(:,l);
            Q(:,l)=Q(:,l)-R(i,l)*Q(:,i);
        end
    end
    Q1=Q(1:NR,1:NT);
    Q2=Q(1:NT,1:NT);
   
    y=Q1'*x(:,j);
    c_u=zeros(NT,1);  %未置换的c矩阵
    %计算第NT个信号的大小
    c_u(NT)=y(NT)/R(NT,NT);
    c_u(NT)=(c_u(NT)>=0)-(c_u(NT)<0)+0; 

    for k=NT-1:-1:1
        d=0;
        for i=k+1:NT
            d=d+R(k,i)*c_u(i);
        end
        z=y(k)-d;
        z=z/R(k,k);
        c_u(k)=(z>=0)-(z<0)+0;  
    end
    
    for h=1:NT
        c(S(h),j)=c_u(h);
    end  
end
c=(c+1)/2;
end


