function e = cel_02(p,u,tm,ym)
    x0 = p(1:3); b = p(4); d = p(5); tf=tm(end);
    [~,x]=rk4_01(x0,u,tf,b,d);
    e = norm(x(:,1)-ym(:,1));
end

