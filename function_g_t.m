function g = function_g_t(x,y,t)
tol = 1e-10;

if abs(t) < tol              % t≈0，区域内部，取初始条件
    g = exp(x + y);
else
    if abs(x - 0) < tol              % 左边界 x = 0
        g = exp(y + t);
    elseif abs(x - 2) < tol          % 右边界 x = 2
        g = exp(2 + y + t);
    elseif abs(y - 0) < tol          % 下边界 y = 0
        g = exp(x + t);
    elseif abs(y - 1) < tol          % 上边界 y = 1
        g = exp(x + 1 + t);
    else
        error('function_g：该点不在边界上且不是初始时刻，本函数不用于内部点 t>0。');
    end
end
end
