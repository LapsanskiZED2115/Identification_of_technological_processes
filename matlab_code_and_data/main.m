close all
clear
%% Charakterystyka statyczna silnika
u_wew = -0.9:0.1:0.9;
c = 0.52;
napiecie = [-1.62,-1.52,-1.4,-1.26,-1.12,-0.98,-0.8,-0.59,-0.34,0,0.35,0.62,0.83,1.01,1.15,1.3,1.43,1.54,1.64];
analog_input = [-3.33,-3.13,-2.85,-2.57,-2.3,-1.98,-1.64,-1.22,-0.69,0,0.7,1.26,1.69,2.05,2.35,2.65,2.9,3.14,3.3];
rpm = napiecie./c.*1000;
wsp = polyfit(rpm,u_wew,3);
analog_wsp = polyfit(analog_input,rpm,1);

figure;
grid on; hold on;
title('Sterowanie obrotami (napiecie)'); ylabel('Sterowanie'); xlabel('Obroty [rpm]')
plot(rpm,u_wew,'o-'); plot(rpm,polyval(wsp,rpm));
legend('Pomiary','Regresja wielomianowa')

figure;
plot(polyval(analog_wsp,analog_input),analog_input)
grid on
title('Sterowanie obrotami (analog input)')
ylabel('Analog input'); xlabel('Obroty [rpm]')
%% Wyznaczanie K dla silnika
load dane

t = u.time;
u = u.signals.values;
y = y.signals.values;

x0 = y(1); K_pocz = 1000;
cel = @(K) cel_00([x0,K],u,t,y);
K_min = fminsearch(cel,K_pocz);

[t,x]=rk4(0,u,t(end),K_min);
figure; plot(t,x,t,y); grid on;title("Wyznaczenie K dla silnika"); legend('Model','Rzeczywiste')
%% Siła ciągu
sila = [-62,-53,-48,-40,-34,-26,-20,-14,-7,-2,0,3,9,18,26,38,50,56,68,80,90,95]*9.81/1000;
rpm = [-3217,-3000,-2815,-2655,-2406,-2170,-1885,-1548,-1135,-655,0,655,1150,1545,1875,2180,2455,2533,2814,2970,3140,3240];

wsp_ciag = polyfit(rpm,sila,3);

figure; grid on; hold on 
plot(rpm,polyval(wsp_ciag,rpm)); plot(rpm,sila);
title('Siła ciągu śmigła'); ylabel('Siła [N]'); xlabel('Obroty [rpm]')
legend('Model','Rzeczywiste',Location='northwest')
%% Równanie wychylenia dla zerowego sterowania
load angle_helikopter.mat
angle = angle_helikopter.signals(1).values;
angle = angle(100:end);
tm = angle_helikopter.time;
tm = tm(100:end);

x0 = [0;0]; T0 = tm(2)-tm(1);
x0(1) = angle(1); x0(2) = (angle(2)-angle(1))/(T0);
a=0.01;c=5;al=0.49;
LB = [-100 -100 0 3 0]'; UB = [100 100 2 15 15]';
options = optimset('display','iter');
[xopt,resnorm,resiudal,exitflag,output,lambda,jacobian] =...
lsqnonlin('cel_01',[x0;a;c;al],LB,UB,options,tm,angle);
x0=xopt(1:2);a=xopt(3); c=xopt(4); al=xopt(5); 

[t,x] = ode45(@(t,x) rhs_01(t,x,a,c,al),tm,x0);
figure; plot(t,x(:,1),tm,angle); grid on; legend("Model","Rzeczywiste")
%% Parametr kp
load oscylacje_ze_sterowaniem.mat
u = u_signal.signals(1).values;
ym = y_signal.signals(1).values;
tm = y_signal.time;
x0 = [0;0;0];

k=12; ym = ym(k:end,:); tm = tm(k:end,1); u = u(k:end,1);
T0=tm(2)-tm(1);

n0=1; n1=length(tm); tm=tm(n0:n1); ym=ym(n0:n1,:); tm=tm-tm(1); u=u(n0:n1);
x0(1,1)=ym(1,1); x0(2,1)=(ym(5,1)-ym(1,1))/(5*T0); x0(3,1)=ym(1,2);

b=0*0.0107388; d=4.5646;
LB=[-100 -1000 -20000 -10 -30]';UB=[100 1000 20000 10 500]';
options=optimset('display','iter');

[xopt,resnorm,residual,exitflag,output,lambda,jacobian]=...
lsqnonlin('cel_02',[x0;b;d],LB,UB,options,u,tm,ym);

x0=xopt(1:3);b=xopt(4);d=xopt(5);
[t,x]=rk4_01(x0,u,tm(end),b,d);
t = (0:0.01:t)';
figure;plot(t,x(:,1),tm,ym(:,1));grid on; legend("Model","Rzeczywiste")
%% Regulator LQ

[xr,ur] = get_eqpoint(0);
[A,B,C,D] = get_abcd(xr);
G=[0 0;1 0;0 1];
sg1=1e-1;sg2=0.5;Q=diag([sg1^2 sg2^2]);
sgv1=1e-1;sgv2=0.02;R=diag([sgv1^2 sgv2^2]);

% LQ controller 
L=lqe(A,G,C,Q,R);
Q=diag([1 100 1]);R=1;
K=lqr(A,B,Q,R);% LQ gain
% LQ controller with integrator
A1=[A zeros(3,1);1 zeros(1,3)];B1=[B;0]; 
Q1=diag([100000 100 0.01 100]);R1=1e4;
K=lqr(A1,B1,Q1,R1);% LQ with integrator gain
%%
load 'pomiary_lq (1).mat'

x3_pomiar = x3_pomiar1.signals.values;
x3_kalman = x3_kalman2.signals.values;
calka_bledu_ = calka_bledu.signals.values;

x2_kalman = x2_kalman1.signals.values;
u_lq = u_lq.signals.values;
t = x1_pomiar.time;

x1_pomiar_val = x1_pomiar.signals.values;
x1_kal = x1_kalman.signals.values;

% Subploty dla x1, x2 i x3
figure;

% Wykres dla x1 (Położenie śmigła)
subplot(3, 1, 1);
plot(t, x1_kal, 'LineWidth', 2); % Estymowana wartość
hold on;
plot(t, x1_pomiar_val, 'LineWidth', 2); % Rzeczywisty pomiar
grid on;
title("Położenie śmigła x1 ");
xlabel("Czas [s]");
ylabel("wychylenie [°]");
legend("Estymowana wartość", "Rzeczywisty pomiar");

% Wykres dla x2 (Prędkość belki)
subplot(3, 1, 2);
plot(t, x2_kalman, 'LineWidth', 2); % Estymowana wartość
grid on;
title("Prędkość belki");
xlabel("Czas [s]");
ylabel("ω[rad/s]");
legend("Estymowana wartość");

% Wykres dla x3 (Prędkość śmigła)
subplot(3, 1, 3);
plot(t, x3_kalman, 'LineWidth', 3); % Estymowana wartość
hold on;
plot(t, x3_pomiar, 'LineWidth', 2); % Rzeczywisty pomiar
grid on;
title("Prędkość śmigła");
xlabel("Czas [s]");
ylabel("ω[RPM]");
legend("Estymowana wartość", "Rzeczywisty pomiar");
hold off;

% Wykres sterowania
figure;
plot(t, u_lq, 'LineWidth', 2);
grid on;
title("Sterowanie");
xlabel("Czas [s]");
ylabel("u");
legend("Sterowanie");

% Wykres całki błędu
figure;
plot(t, calka_bledu_, 'LineWidth', 2);
grid on;
title("Całka błędu");
xlabel("Czas [s]");
ylabel("Całka błędu");

