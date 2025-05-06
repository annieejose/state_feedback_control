clear all 
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
%Euler integration 
%Initializing variables 
t0 = 0; 
tf = 10; 
time = t0:0.01:tf;
X = [30;0]; 
T = t0; 
X_ARR = zeros(length(time),2); % Preallocate the solution matrix 
X_ARR(1,:) = X; % Set the initial condition 
T_ARR = zeros(length(time),1); % Preallocate the solution matrix 
T_ARR(1,1) = t0; % Set the initial condition 
K = [1.5999 2.39506]; 
%Euler integration 
for k=2:1:length(time) 
[X,T,U] = eulr1(X,K,T,M,N,B); 
X_ARR(k,:) = reshape(X,1,2); 
T_ARR(k,:) = T; 
U_ARR(k-1,:) = U; 
end 

figure(3);
plot(T_ARR,X_ARR);hold on;
title('States trajectory');
ylabel('X');
xlabel('time');
legend('X_1','X_2');

%Current
Current = U_ARR ;
figure(5);
plot(T_ARR(1:length(time)-1),Current);
title('Control input trajectory');
hold on;
ylabel('Control input');
xlabel('time');

function [ x,t,u ] = eulr1(x0,P,t0,M,N,B) 
h=0.01; 
[k,u]=state_equation1(x0,P,M,N,B); 
x = x0  + k.*h; 
t = t0 + h; 
end 
function [XDOT,U] = state_equation1(X,K,M,N,B) 
U = N*sin(X(1)) + B*X(2) + M*(-K(1)*X(1)-K(2)*X(2)); 
X1DOT =  X(2); 
X2DOT = (1/M)*(U - N*sin(X(1)-B*X(2))); 
XDOT = [X1DOT; X2DOT]; 
end