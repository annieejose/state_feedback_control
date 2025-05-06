clear all 
close all
clc 
A1 = [0,1 ,0
-14.3634920124791,-0.0208897982795376, 0 
1 0 0]; 
B1 = [0 
1.15697344317439 
0]; 
C1 = [1 0 0]; 
D1 = 0; 
mos = 0.15; 
t_s = 5; 
Seta = sqrt((log(mos))^2/(pi^2+(log(mos))^2)) 
Wn = 4/(Seta*t_s) 
s1 = -Seta*Wn+1i*sqrt(1-Seta^2)*Wn; 
s2 = -Seta*Wn-1i*sqrt(1-Seta^2)*Wn; 
p2 = [s1 s2 -30]; 
K1 = place(A1,B1,p2);%Controller gain 
tf4 = ss(A1,B1,C1,D1); 
figure(5); 
step(tf4); 
A_a1 = A1-B1*K1; 
tf3 = ss(A_a1, B1, C1, D1); 
figure(1); 
step(tf3); 
t0 = 0; 
tf = 10; 
time = t0:0.01:tf; 
X = [30;0;0]; 
T = t0; 
X_ARR = zeros(length(time),3); % Preallocate the solution matrix 
X_ARR(1,:) = X; % Set the initial condition 
T_ARR = zeros(length(time),1); % Preallocate the time matrix 
T_ARR(1,1) = t0; % Set the initial condition 
for k=2:1:length(time) 
[X,T,U] = eulr(X,K1,T,A1,B1); 
X_ARR(k,:) = reshape(X,1,3); 
T_ARR(k,:) = T; 
U_ARR(k-1,:) = U; 
end 
figure(3); 
plot(T_ARR,X_ARR);hold on; 
title('State trajectories'); 
ylabel('X'); 
xlabel('time'); 
legend('X_1','X_2','X_3'); 
% Current 
Current = U_ARR ; 
figure(6); 
plot(T_ARR(1:length(time)-1),Current); 
title('Control input history'); 
hold on; 
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

