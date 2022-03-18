% AMATH 301
% University of Washington
% HW 6
% Due 11/24/2021

clear all; close all; clc

% Define Parameters
M = 200;    % rows
N = 10;     % columns
x = linspace(0,2*pi,M);

%% Part (a)
% Create matrix of orthonormal vectors of sin bases (A1)

A1 = zeros(M,N);

% Debug plot (next 2 rows)
% figure(1);
% hold on;

for n = 1:N
    sn = sin(n*x/2);
    sn = sn/norm(sn);
    A1(:,n) = sn;
%     plot(x,A1(:,n)); % Debug plot
end

% Create correlation matrix (A2) to verify A1 is orthonormal
A2 = A1'*A1;

% Debug A2 plot
% figure(2)
% imagesc(A2); colorbar;

%% Part (b)
% Create matrix of orthonormal vectors of cos bases (A3)

A3 = zeros(M,N);

% Debug plot (next 2 rows)
% figure(3);
% hold on;

for n = 0:N-1
    pn = cos(n*x/2);
    pn = pn/norm(pn);
    A3(:,n+1) = pn;
%     plot(x,A3(:,n+1));
end

% Create correlation matrix (A4) to verify A3 is orthonormal
A4 = A3'*A3;

% Debug A4 plot
% figure(4)
% imagesc(A4); colorbar;

%% Part (c)
% Create correlation matrix (A5) between A1 and A3, verify not orthonormal

A5 = A1'*A3;

% Debug A5 plot
% figure(5)
% imagesc(A5); colorbar;

%% Part (d)
% Compute error of f(x) = cos(x) and function created from A1

f = cos(x);
A6 = zeros(1,N); % E(n)
ansn = 0;
for n = 1:N
    % sn = A1(:,n);
    % an = dot(f,A1(:,n));
    ansn = ansn + dot(f,A1(:,n))*A1(:,n); % sum(an dot sn)
    A6(n) = norm(f' - ansn); % E(n)
end
