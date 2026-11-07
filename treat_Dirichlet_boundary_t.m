function [A,b] = treat_Dirichlet_boundary_t(A, b, boundarynodes, Pb, g_t,t, dirichlet_tag)
% 按“行置 1 法”处理 Dirichlet 边界条件
% g 可以是函数句柄，也可以是函数名字符串，按 g(x,y) 两个坐标参数调用
% dirichlet_tag 为边界节点标记，默认 -1

if nargin < 6 || isempty(dirichlet_tag)
    dirichlet_tag = -1;
end

% 【修改】原代码第11行误写为 g = str2func(g)，变量 g 未定义；应统一为 g_t
if ischar(g_t) || isstring(g_t)
    g_t = str2func(g_t);
end

nbn = size(boundarynodes, 2);

for k = 1 : nbn

    if boundarynodes(1, k) == dirichlet_tag

        i = boundarynodes(2, k);

        A(i, :) = 0;          % 边界节点所在行清零
        A(i, i) = 1;          % 对角元置 1
        b(i) = g_t(Pb(1,i), Pb(2,i),t);

    end
end

end