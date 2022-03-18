% AMATH 301
% University of Washington
% HW 4
% Due 11/4/2021

clear all; close all; clc

% define parameters
dt = 0.01;  % given time step
t = 0:dt:4; % given time interval

% Given equation: v(t) = exp(-0.1t)cos(5t) + t^2 - 0.1t^4
v = exp(-0.1.*t).*cos(5.*t) + t.^2 - 0.1.*t.^4;

% DEBUG PLOTTING %
% 
% % Exact integral
% x_exact = exp(-0.1*t).*(500.*sin(5.*t)-10.*cos(5.*t))/2501 - (t.^5)/50 + (t.^3)/3;
% 
% f1 = figure('Name','Integration Plots');
% 
% plot(t,v); hold on;
% plot(t,x_exact,'LineWidth',[2]);

%% Problem 1a

% A1: integral of v(t) over given t using trapezoidal rule
A1(1) = 0;
for j=2:length(t)
  A=dt*(v(j-1)+v(j))/2;   
  A1(j)=A+A1(j-1);
end

% A2: integral of v(t) over given t using Simpson's rule
A2(1)=0; count=1;
for j=1:2:length(t)-2
  A=dt*(v(j)+4*v(j+1)+v(j+2))/3;   
  A2(count+1)=A+A2(count);
  count=count+1;
end

% DEBUG PLOTTING %
% 
% plot(t,A1,'mo','LineWidth',[1])
% plot(t(1:2:end),A2,'go','LineWidth',[1])

%% Problem 1b

% calc derivative w/ 2nd order center-differencing scheme.
% Use forward- and backward-differencing for endpoints.

A3(1) = (-3*v(1)+4*v(2)-v(3))/(2*dt); % fwd-difference for first point

for j=2:length(t)-1
   A3(j)=( v(j+1)-v(j-1) )/(2*dt);
end

% backward-difference for last point
A3(length(t)) = (3*v(length(t))-4*v(length(t)-1)+v(length(t)-2))/(2*dt);

% DEBUG PLOTTING %
% 
% % Exact derivative of v
% a_exact = -5.*exp(-0.1.*t).*sin(5.*t)-0.1.*exp(-0.1.*t).*cos(5.*t) ...
%     -(2.*t.^3)/5 + 2.*t;
% 
% f2 = figure('Name','2nd Order diff plots');
% plot(t,v); hold on;
% plot(t,a_exact,'g','LineWidth',[2]);
% plot(t,A3,'mo','LineWidth',[1]);

%% Problem 1c

% Calc derivative w/ 4th order center-diff 
% Use 2nd order f- b-diff for 4 end points

A4(1) = (-3*v(1)+4*v(2)-v(3))/(2*dt); % fwd-difference for 1st point
A4(2) = (-3*v(2)+4*v(3)-v(4))/(2*dt); % fwd-difference for 2nd point

for j=3:length(t)-2
   A4(j)=( -v(j+2) + 8*v(j+1) - 8*v(j-1) + v(j-2) )/(12*dt);
end

% bkwd-difference for next-to-last point
A4(length(t)-1) = (3*v(length(t)-1)-4*v(length(t)-2)+v(length(t)-3))/(2*dt);
% bkwd-difference for last point
A4(length(t)) = (3*v(length(t))-4*v(length(t)-1)+v(length(t)-2))/(2*dt);

% DEBUG PLOTTING %
% 
% f3 = figure('Name','4th order diff plots');
% plot(t,v); hold on;
% plot(t,a_exact,'LineWidth',[2]);
% plot(t,A4,'ro','LineWidth',[1]);

%% Problem 1d

% calculate jerk from A4 with 4th-order accurate center-diff scheme
% calc 4 end points with 2nd-order fwd- bkwd-diff schemes

A5(1) = (-3*A4(1)+4*A4(2)-A4(3))/(2*dt); % fwd-difference for 1st point
A5(2) = (-3*A4(2)+4*A4(3)-A4(4))/(2*dt); % fwd-difference for 2nd point

for j=3:length(t)-2
   A5(j)=( -A4(j+2) + 8*A4(j+1) - 8*A4(j-1) + A4(j-2) )/(12*dt);
end

% bkwd-difference for next-to-last point
A5(length(t)-1) = (3*A4(length(t)-1)-4*A4(length(t)-2)+A4(length(t)-3))/(2*dt);
% bkwd-difference for last point
A5(length(t)) = (3*A4(length(t))-4*A4(length(t)-1)+A4(length(t)-2))/(2*dt);


% DEBUG PLOTTING %
% 
% % Exact 2nd derivative (jerk)
% j_exact = 0.01.*exp(-0.1.*t).*(100.*sin(5.*t)-2499.*cos(5.*t)+...
%     (200-100.*t.^2).*exp(.1.*t));
% 
% f4 = figure('Name','4th order jerk plots');
% plot(t,v); hold on;
% plot(t,j_exact,'LineWidth',[2]);
% plot(t,A5,'co','LineWidth',[1]);

%% Problem 1e

% calculate jerk directly from v with 4th-order accurate center-diff scheme
% calc 4 end points with 2nd-order fwd- bkwd-diff schemes

A6(1) = (2*v(1)-5*v(2)+4*v(3)-v(4))/(dt^2); % fwd-diff for 1st point
A6(2) = (2*v(2)-5*v(3)+4*v(4)-v(5))/(dt^2); % fwd-diff for 2nd point

for j=3:length(t)-2
   A6(j)=(-v(j+2)+16*v(j+1)-30*v(j)+16*v(j-1)-v(j-2))/(12*dt^2);
end

% forward-difference for next-to-last point
A6(length(t)-1) = (2*v(length(t)-1)-5*v(length(t)-2)+4*v(length(t)-3)...
    -v(length(t)-4))/(dt^2);
% forward-difference for last point
A6(length(t)) = (2*v(length(t))-5*v(length(t)-1)+4*v(length(t)-2)...
    -v(length(t)-3))/(dt^2);

% DEBUG PLOTTING %
% 
% plot(t,A6,'go','LineWidth',[1]);

%% Problem 1f

A7 = norm(A5-A6);