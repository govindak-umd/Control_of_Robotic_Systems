%% Part D(b) - Non-Linear Calculations
%clearing all previous outputs
clear all
clc
%Declaring the new output variables
y0 = [5; 0; 30; 0; 60; 0]
tspan = 0:0.01:5000;%defining the timespan
[t1,y1] = ode45(@doublepend,tspan,y0); %using ode45 function, specifically built
%for non-linear systems
plot(t1,y1)%plotting the function output on a 2D graph
grid on% making the grid lines visible

%% Defining doublepend function
function dydt = doublepend(t,y)
M=1000;
m1=100;
m2=100;
l1=20;
l2=10;
g=9.81;
A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];
B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];
Q=[1000 0 0 0 0 0;
   0 1000 0 0 0 0;
   0 0 100 0 0 0;
   0 0 0 100 0 0;
   0 0 0 0 100 0;
   0 0 0 0 0 100];
R=1;
[K_val, P_mat, Poles] = lqr(A,B,Q,R);
F=-K_val*y;
dydt=zeros(6,1);
% y(1)=x; y(2)=xdot; y(3)=theta1;   y(4)=theta1dot;  y(5)=theta2;    y(6)=theta2dot;
dydt(1) = y(2); %XD
dydt(2)=(F-(g/2)*(m1*sind(2*y(3))+m2*sind(2*y(5)))-(m1*l1*(y(4)^2)*sind(y(3)))-(m2*l2*(y(6)^2)*sind(y(5))))/(M+m1*((sind(y(3)))^2)+m2*((sind(y(5)))^2));%xDD
dydt(3)= y(4); %theta 1D
dydt(4)= (dydt(2)*cosd(y(3))-g*(sind(y(3))))/l1'; %theta 1 Ddot;
dydt(5)= y(6); %theta 2D
dydt(6)= (dydt(2)*cosd(y(5))-g*(sind(y(5))))/l2; %theta 2Ddot;
end