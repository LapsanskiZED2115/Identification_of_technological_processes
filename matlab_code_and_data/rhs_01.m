function dx=rhs_01(t,x,a,c,al)
    dx1 = x(2);
    dx2 = -a*x(2)-c*sin(x(1)+al);
    dx = [dx1;dx2];
end