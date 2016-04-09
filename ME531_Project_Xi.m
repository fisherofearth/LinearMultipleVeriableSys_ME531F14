%% Define the equations of mothion and linearization
syms m g aSys  c 
syms d dD dDD th thD thDD
syms tau
syms d0 dD0 dD0 th0 thD0 thDD0

%% Linearization position
d0 = 0.15;
dD0 = 0;
dD0 = 0;
th0 = pi/6;
thD0 = 0;
thDD0 =0;
%% constant parameters
m = 1;
g=9.8;
c = 1;
aSys = 0.5;

%% --- x0  xTarget ---
x0 = [d0; 0; th0; 0];
xTarget = [d0;  0; atan(aSys/g); 0];

%% construct State Space
F1 = m*g*sin(th);
F2 = m*aSys*cos(th);
F3 = -m*g*cos(th)*d;
F4 = m*(d^2) *thDD;
F5 = m*aSys*sin(th)*d;

F1_L = simplify(taylor(F1, th, th0, 'Order',2));
F2_L = simplify(taylor(F2, th, th0, 'Order',2));
F3_L = simplify(taylor(F3, [th d], [th0 d0], 'Order',2));
F4_L = simplify(taylor(F4, [thDD d], [thDD0 d0], 'Order',2));
F5_L = simplify(taylor(F5, [th d], [th0 d0], 'Order',2));

C2 = g*sin(th0) - aSys*cos(th0) - g*cos(th0)*th0 - aSys*sin(th0) * th0;
C4 =(aSys*m/((d0^2)*m))* d0*th0*cos(th0) + g*m*d0*th0*sin(th0)/((d0^2)*m)  +(d0*m/((d0^2)*m))* 2*d0*thDD0;
Ad = (-(d0*m/((d0^2)*m))*(2*thDD0)- (aSys*m/((d0^2)*m))*sin(th0) - (g*m*cos(th0)/((d0^2)*m))); %d
Ath = ((g*m*d0*sin(th0)/((d0^2)*m))- (aSys*m/((d0^2)*m))*d0*cos(th0)) ;%th
AthD = - (c/((d0^2)*m)) ; %thD  

% ------ A B C D U---------
A = [0 1 0 0;
    0 0  g*cos(th0) + aSys*sin(th0) 0;
    0 0 0 1;
    Ad 0 Ath AthD];

B = [0 0;
    C2 0;
    0 0; 
    C4 (1/((d0^2)*m))];
U = [1; tau];
 C = eye(4);
%  C = [
%      1 0 0 0;
%      0 0 0 0;
%      0 0 0 0;
%      0 0 0 0;
%      ];
% C = [1 0 0 0];
D = [0 0 ;0 0;0 0;0 0];
% D = [0 0];
sys = ss(A,B,C,D);

%% Controllability and Observability 
PP = ctrb(sys);
QQ = obsv(sys);
rank(PP)
rank(QQ)

%% ----------- K -------------
real = -10;
img = 5;
p1 = real - img*i;
p2 = real + img*i;
p3 = real;
p4 = real;
K = place(A,B,[p1 p2 p3 p4])

p1 = -20 + 10i;
p2 = -20 - 10i;
p3 = -100;
p4 = -100;
K = place(A,B,[p1 p2 p3 p4])

%% figure
figure
data1 = simout.data(:,1);
data2 = simout.data(:,3);
time  = simout.time;

[hAx,hLine1,hLine2] = plotyy(time, data1, time, data2);
ylabel(hAx(1),'d');
ylabel(hAx(2),'theta');
xlabel('time');
title(['system response' ' (dTarget =' num2str(d0) ', theta =' num2str(th0)  ', aSys=' num2str(aSys)  ', m=' num2str(m) ')']);

