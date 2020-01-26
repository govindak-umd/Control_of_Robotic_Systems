%% Part C - Controllability check

%clearing all the previous outputs
clc
clear
%declaring the symbolic variables
syms M m1 m2 l1 l2 g;
% Creating my linearised state space equation using A and B matrices
A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];
%declaring the B matrix
B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)]
% Writing the controllability matrix as follows
disp("Controllability matrix =");
Ct= [B A*B A*A*B A*A*A*B A*A*A*A*B A*A*A*A*A*B]
disp("Finding the determinant of the controllability matrix=");
%simplifying to check for the condition in which the matrix isn't
%controllable
disp(simplify(det(Ct)));
disp("Displaying the rank of the controllability matrix =");
%calculating the rank of the matrix. The system is controllable only if the
%controllability matrix is full rank, that is RANK:6 in this case
rank(Ct)
% From this, we can see that the matrix is invertible for all cases,
% other than the case l1=l2, or for extremely high values of M, m1 
% m2, which is not a realistic case. Therefore, we can conclude that 
% for l1=l2, the matrix is not controllable. This can be shown as follows:
disp("for l1 = l2, Controllability matrix is")
Ct1 = subs(Ct,l1,l2) %using the subs function to make l1 = l2
disp("Displaying rank of the new matrix =")
%Calculating the rank in this case. 
rank(Ct1)
%if loop to display the condition of the system in either case
if (rank(Ct1)==rank(Ct))
    disp("System is controllable as new rank and old rank are same")
else
    disp("system is not controllable as new rank and old rank are different")
end