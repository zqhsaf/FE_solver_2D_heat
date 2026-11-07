function [Linf_err, L2_err, H1_err] = compute_error_2D(solution, P, Pb, T, Tb, basis_type, t)
% 【重写】计算 t 时刻有限元数值解相对解析解 u=e^(x+y+t) 的三类误差
%
% 输入：
%   solution   : 有限元解向量（节点值），维度 N×1
%   P, T       : 顶点坐标与单元连接表（用于单元积分）
%   Pb, Tb     : 有限元节点坐标与局部-全局对应表（用于节点误差）
%   basis_type : 基类型（201 一次 / 202 二次）
%   t          : 当前时间（解析解 u=e^(x+y+t) 依赖于 t）
%
% 输出：
%   Linf_err : ||u-u_h||_∞  全部有限元节点上的最大绝对误差
%   L2_err   : ||u-u_h||_0  L2 范数误差（单元高斯积分）
%   H1_err   : |u-u_h|_1    H1 半范数误差（梯度差的 L2 范数）

switch basis_type
    case 201
        number_of_local_basis = 3;
    case 202
        number_of_local_basis = 6;
    otherwise
        error('compute_error_2D：不支持的 basis_type（201 一次，202 二次）。');
end

Gpn = 5;                       % 7 点 5 阶三角高斯公式
number_of_elements = size(T, 2);

L2_num = 0;   % L2 误差分子
H1_num = 0;   % H1 半范数误差分子

for n = 1 : number_of_elements

    vertices = P(:, T(:, n));
    [Gauss_weights, Gauss_nodes] = generate_Gauss_local_tri(vertices, Gpn);
    Ng = length(Gauss_weights);

    for k = 1 : Ng
        x = Gauss_nodes(k, 1);
        y = Gauss_nodes(k, 2);

        % --- 数值解 u_h 及其梯度（局部基线性组合） ---
        uh  = 0;  uhx = 0;  uhy = 0;
        for beta = 1 : number_of_local_basis
            coef = solution(Tb(beta, n));
            uh  = uh  + coef * FE_basis_local_fun_2D(x, y, vertices, basis_type, beta, 0, 0);
            uhx = uhx + coef * FE_basis_local_fun_2D(x, y, vertices, basis_type, beta, 1, 0);
            uhy = uhy + coef * FE_basis_local_fun_2D(x, y, vertices, basis_type, beta, 0, 1);
        end

        % --- 解析解 u 及其梯度（代入当前时间 t） ---
        u  = exact_sol(x, y, t);
        ux = exact_dx_sol(x, y, t);
        uy = exact_dy_sol(x, y, t);

        % --- 累加 L2 与 H1 半范数的平方 ---
        L2_num = L2_num + Gauss_weights(k) * (uh - u)^2;
        H1_num = H1_num + Gauss_weights(k) * ((uhx - ux)^2 + (uhy - uy)^2);
    end
end

L2_err = sqrt(L2_num);
H1_err = sqrt(H1_num);

% --- L∞ 误差：在全部有限元节点 Pb 上比较 ---
exact_nodes = exact_sol(Pb(1, :), Pb(2, :), t).';
Linf_err = max(abs(solution - exact_nodes));

end
