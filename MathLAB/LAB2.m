%Cm2014, Math, LAB2
%Yihan Xiao, Yizhi Dong

clc;clear;close all;

L = 10;
a = 1;
b_para = 3;
Q_0 = 50;
kappa = 0.5;
rho = 1;
C = 1;
T_out = 300;
T_0 = 400;
v_list = [0, 0.1, 0.5, 1];
N_list = [9, 19, 39, 79];


%% v = 0; N = 9, 19, 39, 79
figure
for i = 1:4
    v = v_list(1);
    N = N_list(i);

    a_general = v*rho*C/kappa;
    h = L/(N+1);
    z = [h:h:L-h]';

    A = zeros(N,N);
    b_mat = zeros(N,1);
    for k = 2:N
        A(k,k) = 2;
        A(k,k-1) = -1-a_general*h/2;
        A(k-1,k) = -1+a_general*h/2;
    end
    A(1,1) = 2;

    b_mat = (h^2/kappa)*Q_0*sin((z-a)/(b_para-a)*pi);
    b_mat(z<=a) = 0; b_mat(z>b_para) = 0;
    b_mat(1) = b_mat(1)+T_0+h*a_general*T_0/2;
    b_mat(N) = b_mat(N)+T_out-h*a_general*T_out/2;

    T = A\b_mat;

    Tplot = [T_0;T;T_out];
    zplot = [0:h:L]';
    plot(zplot, Tplot)
    xlabel('z'),ylabel('T')
    hold on
end
legend(["N = 9", "N = 19", "N = 39", "N = 79"])
title("T(z) for v = 0 and different N")

%% N = 79, v = 0, 0.1, 0.5, 1
figure
for i = 1:4
    v = v_list(i);
    N = N_list(4);

    a_general = v*rho*C/kappa;
    h = L/(N+1);
    z = [h:h:L-h]';

    A = zeros(N,N);
    b_mat = zeros(N,1);
    for k = 2:N
        A(k,k) = 2;
        A(k,k-1) = -1-a_general*h/2;
        A(k-1,k) = -1+a_general*h/2;
    end
    A(1,1) = 2;

    b_mat = (h^2/kappa)*Q_0*sin((z-a)/(b_para-a)*pi);
    b_mat(z<=a) = 0; b_mat(z>b_para) = 0;
    b_mat(1) = b_mat(1)+T_0+h*a_general*T_0/2;
    b_mat(N) = b_mat(N)+T_out-h*a_general*T_out/2;

    T = A\b_mat;

    Tplot = [T_0;T;T_out];
    zplot = [0:h:L]';
    subplot(2,2,i)
    plot(zplot, Tplot)
    xlabel('z'),ylabel('T')
    title("T(z) for N = 79, v = " + num2str(v));
end