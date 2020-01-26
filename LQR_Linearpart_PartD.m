%% Part D(a) - LQR design for linear system
%clearing all the previous outputs
clc
clear all

% Given the constraints for the masses and lengths of the strings
M=1000;%Mass of the cart
m1=100;%mass of Pendulum 1
m2=100;%mass of Pendulum 2
l1=20;%length of the string of Pendulum 1
l2=10;%length of the string of Pendulum 2
g=9.81; %declaring the value of the accelertaion due to gravity in m/s^2
% Porting the A and B matrices here
A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];
B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];
% Checking for the controllability of the given system
if (rank(ctrb(A,B))==6)
    disp("Rank of ctrb matches order of A, system is controllable")
else
    disp("Rank of ctrb doesnt matche order of A, system is uncontrollable")
end
% Now, we give the values of the initial condition for our state variables.
% We are using the state variables in the following format: 
% x(t) = [x x_dot theta_1 theta_1_dot theta_2 theta_2_dot]'
% therefore, the initial conditions are as follows.
x_initial = [0;0;30;0;60;0]; 
%initial disp=0, angles for the two pendulums are as described above. 
% We assume the values of Q and R. For Q, we penalize theta's more than x.
Q=[100 0 0 0 0 0;
   0 100 0 0 0 0;
   0 0 1000 0 0 0;
   0 0 0 1000 0 0;
   0 0 0 0 1000 0;
   0 0 0 0 0 1000];
R=0.01;
% The above values of Q and R are a trade-off and we use them both to
% develop a system as per our priorities
% Let us assume that C matrix is a direct representation of the output
% matrix, which makes D=0
C = eye(6);% To form a 6 X 6 identity matrix
D = 0; % Initialising the D matrix to be Zero
% Lets observe the response of the system for the given initial conditions
% for the given initial values
disp("For the given initial values, the state response is as follows:")
sys1 = ss(A,B,C,D);%MATLAB code for calculating the state space representation of the system
figure
initial(sys1,x_initial)%MATLAB inbuilt function to check the initial response of the system
grid on   %grid lines visible

disp("Now, seeing the results using an LQR controller")
[K_val, P_mat, Poles] = lqr(A,B,Q,R);%In-built MATLAB code for LQR Controllers
K_val %computes the K matrix and displays
P_mat %positive definite matrix calculated for the same
Poles %To see the poles of the given equation
sys2 = ss(A-(B*K_val),B,C,D); %Using the K matrix to define ss
figure
initial(sys2,x_initial)
grid on
% We observe from this, for higher values of Q components, the time it
% takes for the system is die out is much lesser than in other cases. Also,
% lower the value of R, the faster the system stabilizes.