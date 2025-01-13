function dx=predkosc(t,x,u)
% Prawa strona równań stanu
dx_1=x(2);
dx_2=-2*0.071788640050882*2.294877485104032-(2.294877485104032.^2)*sin(x(1));
dx = [dx_1;dx_2];