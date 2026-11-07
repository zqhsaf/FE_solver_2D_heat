function uy = exact_dy_sol(x,y,t)
% 解析解 u=e^(x+y+t) 对 y 的偏导数：u_y = e^(x+y+t)
uy = exp(x + y + t);
end
