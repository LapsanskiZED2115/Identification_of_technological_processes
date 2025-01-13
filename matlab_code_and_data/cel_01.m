function e = cel_01(p,tm,ym)
    x0 = p(1:2); a=p(3); c=p(4); al=p(5); tf = tm(end);
    [t,x] = ode45(@(t,x) rhs_01(t,x,a,c,al),tm,x0);
    e = sum((x(:,1)-ym).^2);
end