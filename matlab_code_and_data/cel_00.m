function e = cel_00(p,u,tm,ym)
    x0 = p(1); K = p(2);
    [t,x]=rk4(x0,u,9.9,K);
    e = sum((x(:,1)-ym).^2);
end