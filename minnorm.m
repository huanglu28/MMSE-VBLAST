function [wk,kk] = minnorm(G,MM)

%本函数的作用是寻找G的范数最小的下标数属于集合M的行,M等于MM消除0元素的集合。
%首先处理各种输入情况到标准情况：G是函数作用矩阵，M是可供比较范数的行序号集合,M是输入MM清除0元素的结果
if nargin==1             %是用来判断输入变量个数的函数
    M=[1:size(G,1)];     % size(G,1)返回矩阵G的行数，例如共有4行，则M=[1:4]={1,2,3,4}
else 
    M=[];
    for i=1:length(MM)   %
        if MM(i)~=0
            temp=MM(i);
            M=[M,temp];
        end
    end
end
%用GG存放G的所有行范数
GG=[];
for i=1:size(G,1)
    GG=[GG;norm(G(i,:))];
end 
%在M集合中挑选对应最小范数的行的元素赋给kk，其对应的最小范数行赋给wk。
kk=M(1);
wk=G(M(1),:);
for i=M
    if GG(i)<=GG(kk)
        kk=i;
        wk=G(i,:);
    end
end