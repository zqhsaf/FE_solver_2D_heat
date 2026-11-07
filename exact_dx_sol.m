function ux = exact_dx_sol(x,y,t)
% 解析解 u=e^(x+y+t) 对 x 的偏导数：u_x = e^(x+y+t)
ux = exp(x + y + t);
end
