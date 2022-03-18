% AMATH 301
% University of Washington
% HW 5
% Due 11/11/2021

clear all; close all; clc

% define parameters
dt = 0.5;
n = 30;
tspan = 0:dt:n;
y0 = [0.1; -1]; % [y(0); y'(0)] Initial conditions

%% Problem 1a, A1

[t,A1] = ode45(@vdpe1, tspan, y0);

% Debug plot
% figure(1);
% plot(t,A1(:,1),'-o',t,A1(:,2),'-o')
% title('Solution of van der Pol Equation (\epsilon = 0.1) with ODE45');
% xlabel('Time t');
% ylabel('Solution y');
% legend('y_1','y_2')

%% Problem 1a, A2

[t,A2] = ode45(@vdpe2, tspan, y0);

% Debug plot
% figure(2);
% plot(t,A2(:,1),'-o',t,A2(:,2),'-o')
% title('Solution of van der Pol Equation (\epsilon = 1) with ODE45');
% xlabel('Time t');
% ylabel('Solution y');
% legend('y_1','y_2')

%% Problem 1a, A3

[t,A3] = ode45(@vdpe3, tspan, y0);

% Debug plot
% figure(3);
% plot(t,A3(:,1),'-o',t,A3(:,2),'-o')
% title('Solution of van der Pol Equation (\epsilon = 20) with ODE45');
% xlabel('Time t');
% ylabel('Solution y');
% legend('y_1','y_2')

%% Problem 1b, A4

[t1,A41] = ode45(@vdpe1, tspan(1:10/dt+1), y0); % t1: [0,10]
[t2,A42] = ode45(@vdpe2, tspan(10/dt+1:20/dt+1), A41(end,:)');%[10,20]
[t3,A43] = ode45(@vdpe3, tspan(20/dt+1:end), A42(end,:)'); %[20,30]

A4 = [A41; A42(2:end,:); A43(2:end,:)];

% Debug plot
% figure(4);
% plot(t,A4(:,1),'-o',t,A4(:,2),'-o')
% title('Solution of van der Pol Equation (\epsilon = [0.1,1,20]), with ODE45');
% xlabel('Time t');
% ylabel('Solution y');
% legend('y_1','y_2')

%% Van der Pol Function

% Equation: y1'' - e(1-y1^2)y1' + y1 = 0
% e = some constant > 0

% Substitute y1' = y2 to make new first-order system
% y1' = y2
% y2' = e(1-y1^2)y2-y1

% Modified versions of dydt function defined in ode45() help documentation

function dydt = vdpe1(t,y)
    dydt = [y(2); 0.1*(1-y(1)^2)*y(2)-y(1)];
end

function dydt = vdpe2(t,y)
    dydt = [y(2); 1*(1-y(1)^2)*y(2)-y(1)];
end

function dydt = vdpe3(t,y)
    dydt = [y(2); 20*(1-y(1)^2)*y(2)-y(1)];
end