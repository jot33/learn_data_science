% AMATH 301
% University of Washington
% HW 7
% Due 12/2/2021

clear all; close all; clc;

%% Problem 1a

% Define Parameters
x = 1:24;
y = [75, 77, 76, 73, 69, 68, 63, 59, 57, ...
     55, 54, 52, 50, 50, 49, 49, 49, 50, ...
     54, 56, 59, 63, 67, 72];



pcoeff = polyfit(x,y,2);
yp = polyval(pcoeff,x);
n = length(yp);
A1 = sqrt( sum( ( abs(yp-y) ).^2  )/n ); % E2
xp = 1:0.01:24;
A2 = polyval(pcoeff,xp)';

% Debug Plot
% figure(1)
% plot(x,y,'o'); hold on;
% plot(xp,A2);
% plot(x,yp);
% legend('y','A2','yp');

%% Problem 1b

A3 = interp1(x,y,xp)';
A4 = spline(x,y,xp)';

% Debug Plot
% figure(2);
% plot(xp,A3,'c+'); hold on;
% plot(xp,A4,'r*');
% plot(x,y,'ko');

%% Problem 1c

% Initial guess based on data
a = (max(y)-min(y))/2; % amplitude
b = pi/14; % angular frequency
c = 63; % offset
c0 = [a; b; c];

% anonymous function 
f = @(cin) E2(cin,x,y);

coeffs = fminsearch(f,c0);

A5 = f(coeffs); % E2 with coefficients using minimized error

fk = @(cin,xdata) cin(1)*cos(cin(2)*xdata)+cin(3);
yold = fk(c0,x);
A6 = fk(coeffs,xp)';

% Debug plot
% figure(3);
% plot(x,y,'o'); hold on;
% plot(x,yold,'m'); % sinusoidal with initial guess
% plot(xp,A6','g'); % sinusoidal with minimized error
% legend('data points','yold','ynew');

%% Problem 2

% define velocity vector
v = [ 30, 35, 33, 32, 34, 37, 39, 38, 36, 36, 37, 39, 42, 45, 45,...
      41, 40, 39, 42, 44, 47, 49, 50, 49, 46, 48, 50, 53, 55, 54, 53];
  
t = 0:30;

% Initial guess
a2 = 3;
b2 = pi/4;
c2 = 2/3;
d2 = 32;
c0_2 = [a2,b2,c2,d2];

f2 = @(cin) E2_2(cin,t,v);

coeffs_2 = fminsearch(f2,c0_2);

fk2 = @(cin,xdata) cin(1)*cos(cin(2)*xdata)+cin(3)*xdata+cin(4);
vold = fk2(c0_2,t); % velocity function based on initial guess

tp = 0:0.01:30;

A7 = fk2(coeffs_2,tp);

% Debug plots
% figure(4);
% plot(t,v,'o'); hold on;
% plot(t,vold,'m');
% plot(tp,A7,'g');
% legend('data points', 'v_{old}','v_{new}');

%% Custom functions

% Custom function for finding mean squared error
function E = E2(cin,xdata,ydata)
    A = cin(1);
    B = cin(2);
    C = cin(3);
    n = length(ydata);
    ynew = A*cos(B*xdata)+C;
    E = sqrt(sum((abs(ydata - ynew)).^2)/n);
end

function E = E2_2(cin,xdata,ydata)
    A = cin(1);
    B = cin(2);
    C = cin(3);
    D = cin(4);
    n = length(ydata);
    ynew = A*cos(B*xdata)+C*xdata+D;
    E = sqrt(sum((abs(ydata - ynew)).^2)/n);
end