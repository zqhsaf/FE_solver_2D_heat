function boundaryedges=generate_boundaryedges(N1,N2)
% 一次三角元结构化网格的边界边表（4 行）：
%   第 1 行：边界类型标记，-1 表示 Dirichlet
%   第 2 行：该边界边所属的单元编号
%   第 3、4 行：边界边两个端点的全局节点编号
% 按 bottom -> right -> top -> left 顺序排列，共 2(N1+N2) 条边界边。

N=2*(N1+N2);
boundaryedges=zeros(4,N);
col=0;

% ---------- bottom 边（y=bottom，从左到右），局部边 1-2 ----------
for k=1:N1
    col=col+1;
    a=(k-1)*(N2+1)+1;
    b=k*(N2+1)+1;
    elem=2*(k-1)*N2+1;
    boundaryedges(:,col)=[-1; elem; a; b];
end

% ---------- right 边（x=right，从下到上），局部边 2-3 ----------
for k=1:N2
    col=col+1;
    j=k-1;
    a=N1*(N2+1)+j+1;
    b=N1*(N2+1)+j+2;
    elem=2*((N1-1)*N2+j)+2;
    boundaryedges(:,col)=[-1; elem; a; b];
end

% ---------- top 边（y=top，从右到左），局部边 1-3 ----------
for k=1:N1
    col=col+1;
    i=N1-k;
    a=i*(N2+1)+N2+1;
    b=(i+1)*(N2+1)+N2+1;
    elem=2*(i*N2+(N2-1))+2;
    boundaryedges(:,col)=[-1; elem; a; b];
end

% ---------- left 边（x=left，从上到下），局部边 1-3 ----------
for k=1:N2
    col=col+1;
    j=N2-k;
    a=j+1;
    b=j+2;
    elem=2*j+1;
    boundaryedges(:,col)=[-1; elem; a; b];
end

end