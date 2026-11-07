function [solution, Linf_err, L2_err, H1_err] = FE_solver_2D_heat(start,ending,dt,theta,left,right,bottom,top,basis_type_trial,basis_type_test,N1,N2)
% 二维瞬态热传导方程  u_t - ∇·(2∇u) = f  的三角元有限元求解
% 求解区域：[left,right] × [bottom,top]，四条边均为 Dirichlet 边界
% basis_type_trial / basis_type_test：
%   201 一次三角元（每单元 3 个局部基）；202 二次三角元（每单元 6 个局部基）
% N1,N2：x、y 方向上的均匀剖分段数
%
% 输出：
%   solution : t=ending 时刻的数值解向量
%   Linf_err : ||u-u_h||_∞ 节点最大误差
%   L2_err   : ||u-u_h||_0  L2 范数误差
%   H1_err   : |u-u_h|_1    H1 半范数误差

%% 1. 生成网格与有限元节点
[P,T]=generate_PT(left,right,bottom,top,N1,N2);
[Pb,Tb]=generate_PbTb(left,right,bottom,top,N1,N2,basis_type_trial);

switch basis_type_trial
    case 201
        Pb_trial=Pb; Tb_trial=Tb; number_of_local_basis_fun_trial=3;
    case 202
        Pb_trial=Pb; Tb_trial=Tb; number_of_local_basis_fun_trial=6;
    otherwise
        error('主求解器：不支持的试探基类型（201 一次，202 二次）。');
end

switch basis_type_test
    case 201
        Pb_test=Pb; Tb_test=Tb; number_of_local_basis_fun_test=3;
    case 202
        Pb_test=Pb; Tb_test=Tb; number_of_local_basis_fun_test=6;
    otherwise
        error('主求解器：不支持的检验基类型（201 一次，202 二次）。');
end

%% 2. 生成边界节点（基于坐标，一次/二次元通用）
boundarynodes=generate_boundarynodes(Pb_test,left,right,bottom,top);

matrix_size=[size(Pb_test,2),size(Pb_trial,2)];
number_of_elements=size(T,2);

%% 3. 组装质量矩阵与刚度矩阵
% 【修改】质量矩阵 M = ∫ φ_trial·ψ_test dxdy，试探基与检验基均不求导
M=assemble_matrix_2D('function_one',matrix_size,number_of_elements,P,T, ...
    Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test, ...
    basis_type_trial,0,0,basis_type_test,0,0);

A1=assemble_matrix_2D('function_c',matrix_size,number_of_elements,P,T, ...
    Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test, ...
    basis_type_trial,1,0,basis_type_test,1,0);

A2=assemble_matrix_2D('function_c',matrix_size,number_of_elements,P,T, ...
    Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test, ...
    basis_type_trial,0,1,basis_type_test,0,1); 

A=A1+A2;

% 【修改】θ-离散格式：(M/dt + θA)u^{n+1} = (M/dt − (1−θ)A)u^n + θb^{n+1} + (1−θ)b^n
Atilde=M/dt+theta*A;
Afixed=M/dt-(1-theta)*A;
clear M A;
% 【修改】生成初始向量：传入有限元节点坐标 Pb_trial，返回 t=0 时各节点初始值
X_old=generate_initial_vector(Pb_trial);
number_of_time_step=(ending-start)/dt;

vector_size=[size(Pb_test,2),1];
for m=0:number_of_time_step-1

    tm=m*dt;
    tmp1=(m+1)*dt;

    bm=assemble_vector_2D('function_f',tm,vector_size,number_of_elements,P,T, ...
    Tb_test,number_of_local_basis_fun_test,basis_type_test,0,0);
    bmp1=assemble_vector_2D('function_f',tmp1,vector_size,number_of_elements,P,T, ...
    Tb_test,number_of_local_basis_fun_test,basis_type_test,0,0);
    btilde=theta*bmp1+(1-theta)*bm+Afixed*X_old;
%% 4. 处理 Dirichlet 边界条件
    [Atilde,btilde]=treat_Dirichlet_boundary_t(Atilde,btilde,boundarynodes,Pb_test,'function_g_t',tmp1,-1);
%% 5. 求解整体线性方程组
    solution=Atilde\btilde;
    X_old=solution;
end

%% 6. 计算 t=ending 时刻的三类误差（与解析解 u=e^(x+y+t) 对比）
% 【新增】调用误差计算函数，返回 Linf / L2 / H1 三种误差
[Linf_err, L2_err, H1_err] = compute_error_2D( ...
    solution, P, Pb_trial, T, Tb_trial, basis_type_trial, ending);

end