function c= MMSE_SQRD_PSA(H,x,snr)
% 应用了PSA的MMSE_SQRD算法
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
    %对HH qr分解
    R=zeros(NT,NT);
    Q=HH;
    S=1:NT;
    for i=1:NT
        %R的对角线元素等于qi的1范数
        R(i,i)=norm(Q(:,i),1);
        Q(:,i)=Q(:,i)/R(i,i);
        for l=i+1:NT
            R(i,l)=Q(:,i)'*Q(:,l);
            Q(:,l)=Q(:,l)-R(i,l)*Q(:,i);
        end
    end
    Q1=Q(1:NR,1:NT);
    Q2=Q(1:NT,1:NT);
    
   %Post-Sorting-Algorithm
    k_min=NT;
    for i=NT:-1:2
        error=[];
        for l=1:i
            error=[error,Q2(l,1:i)*Q2(l,1:i)'];
        end
        [~,k_i]=min(error);
        k_min=min([k_min,k_i]);
        if(k_i<i)
            %交换Q2的第i和第k_i行
            temp=Q2(i,:);
            Q2(i,:)=Q2(k_i,:);
            Q2(k_i,:)=temp;
            %交换S的第i列和第k_i列
            temp=S(i);
            S(i)=S(k_i);
            S(k_i)=temp;
        end
        if(k_min<i)
            for(n=1:i)
                %norm_a -- Q2(n,k_min:i)的范数
                a=Q2(n,k_min:i);
                norm_a=sqrt(a*a');
                e_n=zeros(1,(i-k_min+1));
                e_n(i-k_min+1)=1;
                u=(a-norm_a*e_n)/norm(a-norm_a*e_n);
                w=u*a'/(a*u');
                %HR -- householder refletor of Q2(n,k_min:i)
                HR=ones(i-k_min+1,i-k_min+1)-(1+w)*u'*u;
                Q2(n,k_min:i)=Q2(n,k_min:i)*HR;
                Q1(:,k_min:i)=Q1(:,k_min:i)*HR;
            end
        end
    end
    
    R=1/sqrt(1/snr)*inv(Q2);
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

