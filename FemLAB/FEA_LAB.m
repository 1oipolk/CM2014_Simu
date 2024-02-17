% CM2014, FEA, LAB
% Yihan Xiao, Yizhi Dong
clc; clear; close all

% Radius
r = 0.05; % m
% Young's Modulus
E = 100e6; % Pa
% Rotational Inertia
I = pi * r^4 / 4; % kg/m^2
% Cross-section area
A = pi * r^2; % m^2
% Length of the rod
L = 2; % m
% Length of each segment
l = L / 4;
% Initial Angels
beta = deg2rad(70); % rad
lambda = deg2rad(55); % rad
delta = deg2rad(30); % rad
% External Force
F = 20; % N

% Local Stiffness matrix
k_ul = blkdiag(E*A/l, [12*E*I/l^3 6*E*I/l^2; 6*E*I/l^2 4*E*I/l]);
k_ll = blkdiag(-E*A/l, [-12*E*I/l^3 6*E*I/l^2; -6*E*I/l^2 2*E*I/l]);
k_ur = blkdiag(-E*A/l, [-12*E*I/l^3 -6*E*I/l^2; 6*E*I/l^2 2*E*I/l]);
k_lr = blkdiag(E*A/l, [12*E*I/l^3 -6*E*I/l^2; -6*E*I/l^2 4*E*I/l]);
k_e = [k_ul k_ll; k_ur k_lr];

K_e = zeros(15);

% Transformation matrix
Angle = [beta, lambda, delta, 0];

for i = 0:3
    % Transformation matrix
    T = axisT(Angle(i+1));
    % Aixs transform
    K = inv(T)*k_e*T;
    % Global stiffness matrix
    K_e = K_e + blkdiag(zeros(i*3), K, zeros((3-i)*3));
end

%
K_e(:,1:3) = 0;
K_e(1:3,:) = 0;

% Global force vector
F_e = [0; 0; zeros(11, 1); -20; 0];

% Global deformation matrix
D_e = K_e\F_e

X_0 = [0]; Z_0 = [0];
x_i = 0; z_i = 0;

for i = 1:4
    x_i = x_i + l*cos(Angle(i));
    z_i = z_i + l*sin(Angle(i));
    X_0 = [X_0, x_i];
    Z_0 = [Z_0, z_i];
end

X_1 = X_0; Z_1 = Z_0;

for i = 1:4
    X_1(i+1) = X_0(i+1) + D_e(3*i+1);
    Z_1(i+1) = Z_0(i+1) + D_e(3*i+2);
end

plot(X_0, Z_0)
hold on
plot(X_1, Z_1)
title("Change of the Beam")
legend(["Oringinal Plot", "With Deformation"], 'Location', 'southeast')

function T_mat = axisT(angle)
    aixsMat = [cos(angle) sin(angle); -sin(angle) cos(angle)];
    T_mat = blkdiag(aixsMat, 1, aixsMat, 1);
end