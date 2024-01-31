%Cm2014, Math, LAB1
%Yihan Xiao, Yizhi Dong

clc;clear;close all;
%% Theoretical Plot
omega=2;
h=0.1;N=100;
t=[0:h:N*h];
x=cos(omega*t);
xp=-omega*sin(omega*t);
E=omega^2*(x.*x)+xp.*xp;
X=[x' xp'];

figure
subplot(3,2,[1,2])
plot(t,x,t,xp)
xlabel('t'),ylabel('x,xp')
title('Harmonic oscillator, exact trajetories')

subplot(3,2,3)
plot(x,xp)
axis('equal')
xlabel('x'),ylabel('dx/dt')
title('Harmonic oscillator, exact phase portrait')
grid

subplot(3,2,[5,6])
plot(t,E)
xlabel('t')
ylabel('E')
title('Harmonic oscillator, exact Energi vs time')

%% 1) Symplectic Euler
y=[1,0]';t=0;
A=[0 1;-omega^2 0];
Y=[y'];T=[t];E=[omega^2];
for k=1:N
    y=(eye(2)+h*A)*y;
    Y=[Y;y'];
    t=t+h;
    T=[T;t];
    E=[E;omega^2*y(1)*y(1)+y(2)*y(2)];
end
Y_SE=Y;

figure
subplot(3,2,[1,2])
plot(T,Y,'o')
xlabel('t'),ylabel('y1,y2')
title('Harmonic Oscillator, Symplectic Euler trajectories')

subplot(3,2,3)
plot(Y(:,1),Y(:,2),'o')
axis('equal')
xlabel('y1'),ylabel('y2')
title('Harmonic Oscillator, Symplectic Euler phase portrait')

subplot(3,2,[5,6])
plot(T,E)
xlabel('t'),ylabel('E')
title('Harmonic Oscillator, Symplectic Euler Energy vs time')

%% 2) Implicit Midpoint Method
y=[1,0]';t=0;
Y=[y'];T=[t];E=[omega^2];
for k=1:N
    y=(eye(2)-h*A/2)\((eye(2)+h*A/2)*y);
    Y=[Y;y'];
    t=t+h;
    T=[T;t];
    E=[E;omega^2*y(1)*y(1)+y(2)*y(2)];
end
Y_IMM=Y;

figure
subplot(3,2,[1,2])
plot(T,Y,'o')
xlabel('t'),ylabel('y1,y2')
title('Harmonic Oscillator, Implicit Midpoint trajectories')

subplot(3,2,3)
plot(Y(:,1),Y(:,2),'o')
axis('equal')
xlabel('y1'),ylabel('y2')
title('Harmonic Oscillator, Implicit Midpoint phase portrait')

subplot(3,2,[5,6])
plot(T,E)
xlabel('t'),ylabel('E')
title('Harmonic Oscillator, Implicit Midpoint Energy vs time')

%% 3) Verlet’s method
y=[1,0]';t=0;
Y=[y'];T=[t];E=[omega^2];

u1 = 1 + h*0 - omega^2*h^2*1/2;
u2 = 0;

for k=1:N-1
    u1_p1 = (2 - h^2*omega^2)*u1 - y(1);
    u2 = (u1_p1 - y(1))/(2*h);
    y=[u1, u2]';
    Y=[Y;y'];
    t=t+h;
    T=[T;t];
    E=[E;omega^2*y(1)*y(1)+y(2)*y(2)];

    u1 = u1_p1;
end
Y_VM=Y;

figure
subplot(3,2,[1,2])
plot(T,Y,'o')
xlabel('t'),ylabel('y1,y2')
title('Harmonic Oscillator, Verlet’s method trajectories')

subplot(3,2,3)
plot(Y(:,1),Y(:,2),'o')
axis('equal')
xlabel('y1'),ylabel('y2')
title('Harmonic Oscillator, Verlet’s method phase portrait')

subplot(3,2,[5,6])
plot(T,E)
xlabel('t'),ylabel('E')
title('Harmonic Oscillator, Verlet’s method Energy vs time')
