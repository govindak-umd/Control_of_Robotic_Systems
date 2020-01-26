%% Part F - Luenberger observer for linear systems
%clearing all the previous outputs
clc
clear all
%Substituting values for our M, m1, m2, l1 and l2
M=1000;%Mass of the cart
m1=100;%mass of Pendulum 1
m2=100;%mass of Pendulum 2
l1=20;%length of the string of Pendulum 1
l2=10;%length of the string of Pendulum 2
g=9.81; %declaring the value of the accelertaion due to gravity in m/s^2
%Defining our matrices as follows
A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];
B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];
% From previous case, we have determined that only C1, C3 and C4 were
% observable. Hence, we are going to consider only those 3 cases.
C1 = [1 0 0 0 0 0];  %Corresponding to x component
C3 = [1 0 0 0 0 0; 0 0 0 0 1 0]; %cooresponding to x and theat2
C4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]; %cooresponding to x, theta1 and theat2
D = 0;%declarring the D matrix to be zero
% Cosidering the same Q and R matrices chosen before in our code
Q=[100 0 0 0 0 0;
   0 100 0 0 0 0;
   0 0 1000 0 0 0;
   0 0 0 1000 0 0;
   0 0 0 0 1000 0;
   0 0 0 0 0 1000];
R=0.01; %%these are the cost variables from LQR
% Initial Conditions for Leunberger observer - 12 state variables, 
% 6 actual + 6 estimates
x0=[0,0,30,0,60,0,0,0,0,0,0,0];
% state variables order = [x,dx,theta_1,dtheta_1,theta_2,dtheta_2,
% estimates taken in the same order]
% For pole placement, lets choose eigen values with negative real part
poles=[-1;-2;-3;-4;-5;-6];
% Calling LQR function to obtain K matrix
K=lqr(A,B,Q,R);
% Framing L for all three cases where output is observable
% Using the pole placement funciton built into MATLAB
L1 = place(A',C1',poles)' %L1 should be a 6x1 matrix
L3 = place(A',C3',poles)' %L3 should be a 6x2 matrix
L4 = place(A',C4',poles)' %L4 should be a 6x3 matrix

Ac1 = [(A-B*K) B*K; % Luenberger A matrix
        zeros(size(A)) (A-L1*C1)];
Bc = [B;zeros(size(B))];% Luenberger B matrix
Cc1 = [C1 zeros(size(C1))];% Luenberger C matrix

Ac3 = [(A-B*K) B*K;% Luenberger A matrix
        zeros(size(A)) (A-L3*C3)];
Cc3 = [C3 zeros(size(C3))];% Luenberger C matrix

Ac4 = [(A-B*K) B*K;% Luenberger A matrix
        zeros(size(A)) (A-L4*C4)];
Cc4 = [C4 zeros(size(C4))];% Luenberger C matrix

sys1 = ss(Ac1, Bc, Cc1,D);%%MATLAB function to output statespace equations
figure % to launch a new figure WINDOW
initial(sys1,x0)%MATLAB inbuilt function to check the initial response of the system
figure
step(sys1)%Gives the step response output

sys3 = ss(Ac3, Bc, Cc3,D);%MATLAB function to output statespace equations
figure 
initial(sys3,x0)%MATLAB inbuilt function to check the initial response of the system
figure 
step(sys3)%Gives the step response output

sys4 = ss(Ac4, Bc, Cc4, D);%%MATLAB function to output statespace equations
figure 
initial(sys4,x0)%MATLAB inbuilt function to check the initial response of the system
figure 
step(sys4)%Gives the step response output

grid on