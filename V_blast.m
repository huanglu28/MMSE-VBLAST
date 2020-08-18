function c =V_blast( H,x )
% V-blast算法检测多层时空信号
% H -- NR*NT维瑞利信道
% x -- 接收信号
% c -- 解码信号
[NR,NT,L]=size(H);
c=zeros(NT,L);
    for j=1:L
        HH=H(:,:,j);
        S=1:NR;%标记被消去的列
        for i=1:NT
            G=pinv(HH);%G为H的二维广义逆矩阵
            [w,k]=minnorm(G);%k为G的最小行范数所在行
            y=G(k,:)*x(:,j); %生成判决统计量
            temp=S(k);
            %对y进行判决，得到解码信号
            c(temp,j)=1*(y>=0)-1*(y<0)+0;
            %从接受信号中消除c，得到修正后的信号
            x(:,j)=x(:,j)-sqrt(1/NR)*HH(:,k)*c(temp,j);
            HH(:,k)=[];
            S(k)=[];
        end
        end
   c=(c+1)/2;
end






