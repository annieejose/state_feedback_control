clear all 
close all
clc 
% Define the state variables and the differential equations 
syms x1 x2 I 
J = 1.625e-3; 
m = 0.506; 
M = 0.434; 
D = 0.305; 
b = 16.25e-3; 
L = 25e-3; 
R = 5; 
k_t = 0.90; 
del = 1; 
d = 1; 
g = 10; 
M = J/k_t + m*d^2/(3*k_t) + M*d^2/k_t + 2*M*del^2/(5*k_t); 
N = m*d*g/(2*k_t) + M*d*g/k_t; 
B = b/k_t; 
f1 = x2; 
f2 = (1/M)*(I - N*sin(x1)-B*x2); 

% Define an initial guess for the equilibrium point 
eq0 = [0, 0]; 

% Set the maximum number of iterations and the tolerance for convergence 
max_iter = 100; 
tol = 1e-6; 

% Perform fixed-point iteration to find the equilibrium point 
eq = eq0; 
for j = 1:max_iter 
eq_new = myfun(eq,M,I,N,B); 
if norm(eq_new - eq) < tol
    break; 
end 
eq = eq_new; 
end 

% Define the state vector 
x = [x1; x2]; 
% Define the matrix of partial derivatives (Jacobian matrix) 
J1 = jacobian([f1; f2], x); 
J2 = jacobian([f1; f2], I); 
% Evaluate the Jacobian matrix at a specific point 
x0 = [0;0]; 
u0 = 0; 
J_eval1 = double(subs(J1, x, x0)); 
J_eval2 = double(subs(J2, I, 0)); 
A = J_eval1; 
B = J_eval2; 
C = [1 0]; 
D = 0 ; 
tf1 = ss(A, B, C, D); 
figure(1); 
step(tf1);hold on; 
grid on; 
mos = 0.15; 
t_s = 5; 
Seta = sqrt((log(mos))^2/(pi^2+(log(mos))^2)) 
Wn = 4/(Seta*t_s) 
s1 = -Seta*Wn+i*sqrt(1-Seta^2)*Wn; 
s2 = -Seta*Wn-1i*sqrt(1-Seta^2)*Wn; 
p = [s1 s2] 
K = place(A,B,p) 
A_a = A-B*K; 
tf2 = ss(A_a, B, C, D); 
figure(2); 
step(tf2);hold on; 
grid on; 
%Initializing variables 
t0 = 0; 
tf = 10; 
time = t0:0.01:tf; 
X = [30;0]; 
T = t0; 
X_ARR = zeros(length(time),2); % Preallocate the solution matrix 
X_ARR(1,:) = X; % Set the initial condition 
T_ARR = zeros(length(time),1); % Preallocate the time matrix 
T_ARR(1,1) = t0; % Set the initial condition 
%Euler integration to solve for X
for k=2:1:length(time) 
[X,T,U] = eulr(X,K,T,A,B); 
X_ARR(k,:) = reshape(X,1,2); 
T_ARR(k,:) = T; 
U_ARR(k-1,:) = U; 
end 
%Plot of the states - angular motor rotor position and velocity
figure(3); 
plot(T_ARR,X_ARR);hold on;grid on; 
title('States vs t'); 
ylabel('X'); 
xlabel('time'); 
legend('X_1','X_2'); 
%Plot of the control input-Current 
Current = U_ARR ; 
figure(4); 
plot(T_ARR(1:length(time)-1),Current); 
title('Control input history over time'); 
hold on;grid on; 
ylabel('Control input'); 
xlabel('time'); 

function f = myfun(x,M,I,N,B)
I = 0;
f(1) = x(2);
f(2) = (1/M)*(I - N*sin(x(1)-B*x(2)));
end

%Euler integration equation
function [ x,t,u ] = eulr(x0,P,t0,A,B) 
h=0.01; 
[k,u]=state_equation(x0,P,A,B); 
x = x0  + k.*h; 
t = t0 + h; 
end 

%State equation
function [XDOT,U] = state_equation(X,K,A,B) 
U=-K*X; 
XDOT = A*X+B*U; 
end