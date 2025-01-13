function dx=rhs_00(t,x,u,K)
% Prawa strona równań stanu
H= 1.28308710741392e-11*x.^3-9.47671659934249e-10*x.^2+0.000162289932269009*x-0.00418773065630636;%static
dx=K*(u-H);