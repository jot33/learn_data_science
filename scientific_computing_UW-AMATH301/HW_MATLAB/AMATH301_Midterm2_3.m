% AMATH 301
% University of Washington
% Midterm 2
% Due 11/19/2021

clear all; close all; clc

% define parameters
dt = 0.01;
tspan = 0:dt:10;
gamma = 0.1;
y0 = [1; 0.1]; % initial conditions

%% Part a (Euler approximation of pendulum diff eqn)

% Let theta = y1
% f1(tn,yn) = y1' = y2
% f2(tn,yn) = y2' = -gamma*y2 - sin(y1)
% y(t)  = [theta; theta']
%       = [y1_n + dt*y2_n;
%          y2_n + dt*(-gamma*y2_n - sin(y1_n))]
A1 = zeros(2,length(tspan));
A1(:,1) = y0;
for n=1:length(tspan)-1
    A1(:,n+1) = [A1(1,n) + dt*A1(2,n);
                 A1(2,n) + dt*(-gamma*A1(2,n) - sin(A1(1,n)))];
end

%{
ODE45 for comparison to euler method
pend = @(t,y,g) [y(2); -g*y(2)-sin(y(1))];
[t,y45] = ode45(pend, tspan, y0, g=gamma);

% Debug plot
figure(1);
plot(t,y45(:,1),'-o',t,y45(:,2),'-o')
title('Solution of pendulum with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('\theta','\theta''')

% Debug plot
figure(2);
plot(t,A1(1,:)','-o',t,A1(2,:)','-o')
title('Solution of pendulum with Euler');
xlabel('Time t');
ylabel('Solution y');
legend('\theta','\theta''')

% Debug plot
figure(3);
plot(t,A1(1,:)','-o',t,y45(:,1),'-o')
title('Solution of pendulum displacement, Euler vs ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('\theta_{euler}','\theta_{ODE45}')
%}

%% Part b (Trapezoidal rule vector)

A2 = zeros(1,length(tspan));
for n=1:length(tspan)-1
    A2(n+1) = A2(n) + dt*(A1(1,n)+A1(1,n+1))/2;
end

% int_theta = trapz(tspan, A1(1,:)); % for comparison
%% Part b (Simpsons rule vector)

A3 = zeros(1,501);
for j=1:2:length(tspan)-2   
  n = (j+1)/2;
  A3(n+1)= A3(n) + dt*(A1(1,j)+4*A1(1,j+1)+A1(1,j+2))/3;
end

%% Part c (2nd-order accuracy scheme,accel from velocity)
A4 = zeros(1,length(tspan));
A4(1) = (-3*A1(2,1)+4*A1(2,2)-A1(2,3))/(2*dt); % fwd-difference first point

for n=2:length(tspan)-1
   A4(n)=( A1(2,n+1) - A1(2,n-1) )/(2*dt);
end

% backward-difference for last point
A4(1001) = (3*A1(2,1001)-4*A1(2,1001-1)+A1(2,1001-2))/(2*dt);

% debug plot comparison between part a and part c
% plot(tspan,A4,'-o',tspan,A1(2,:),'-o')

%% Part c (4th-order accuracy scheme,accel from velocity)

A5 = zeros(1,length(tspan));
A5(1) = (-3*A1(2,1)+4*A1(2,2)-A1(2,3))/(2*dt); % fwd-difference 1st point
A5(2) = (-3*A1(2,2)+4*A1(2,3)-A1(2,4))/(2*dt); % fwd-difference 2nd point

for n=3:length(tspan)-2
   A5(n)=( -A1(2,n+2) + 8*A1(2,n+1) - 8*A1(2,n-1) + A1(2,n-2) )/(12*dt);
end

% bkwd-difference next-to-last point
A5(1001-1) = (3*A1(2,1001-1)-4*A1(2,1001-2)+A1(2,1001-3))/(2*dt);
% bkwd-difference last point
A5(1001) = (3*A1(2,1001)-4*A1(2,1001-1)+A1(2,1001-2))/(2*dt);

% debug plot comparison between part a and part c
% plot(tspan,A5,'-o',tspan,A1(2,:),'-o')

%% Part d (4th-order accuracy scheme for second derivative)

A6 = zeros(1,length(tspan));
% fwd-difference 1st point
A6(1) = (2*A1(1,1)-5*A1(1,2)+4*A1(1,3)-A1(1,4))/(dt^2);
% fwd-difference 2nd point
A6(2) = (2*A1(1,2)-5*A1(1,3)+4*A1(1,4)-A1(1,5))/(dt^2);

for n=3:length(tspan)-2
   A6(n)=(-A1(1,n+2)+16*A1(1,n+1)-30*A1(1,n)+...
       16*A1(1,n-1)-A1(1,n-2))/(12*dt^2);
end

% bkwd-difference next-to-last point
A6(1001-1) = (2*A1(1,1001-1)-5*A1(1,1001-2)+...
    4*A1(1,1001-3)-A1(1,1001-4))/(dt^2);
% bkwd-difference last point
A6(1001) = (2*A1(1,1001)-5*A1(1,1001-1)+...
    4*A1(1,1001-2)-A1(1,1001-3))/(dt^2);

% debug plot comparing A4,A5,A6 (all methods for accel)
% plot(tspan,A4,'-o',tspan,A5,'-o',tspan,A6,'-o')
% legend('\theta''''_{A4}','\theta''''_{A5}','\theta''''_{A6}')

%% Part e (Adams-Bashforth scheme)

A7 = zeros(2,length(tspan));
A7(:,1:2) = [y0 A1(:,2)]; % use part a for second point

for n=2:length(tspan)-1
    f2n = -gamma*A7(2,n) - sin(A7(1,n)); %f2(t_n,y_n)
    f2n1 = -gamma*A7(2,n-1) - sin(A7(1,n-1)); %f2(t_{n-1},y_{n-1})
    A7(:,n+1) = [A7(1,n) + dt*(3*A7(2,n)-A7(2,n-1))/2;
                 A7(2,n) + dt*(3*f2n - f2n1)/2];
end

% debug plot comparison between A1 and A7
% plot(tspan,A1(1,:),'-o',tspan,A7(1,:),'-o')
% legend('\theta_{A1}','\theta_{A7}')

%% Part f (2nd-order Runge-Kutta scheme using ode23)

pend = @(t,y) [y(2); -gamma*y(2)-sin(y(1))];
[t,y23] = ode23(pend, tspan, y0);
A8 = y23';

% debug plot comparison between A1 and A8
% plot(tspan,A1,'-o',tspan,A8,'-o')
% legend('\theta_{A1}','\theta_{A1}','\theta_{A1}','\theta_{A8}')