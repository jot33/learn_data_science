% AMATH 301
% Midterm 1
% University of Washington
% Due 10/30/2021

clear all; close all; clc

% Debug display (compare for convergence)
% A = [8 2 1; -1 5 1;  1 1 4];
% B = [A(3,:); A(2,:); A(1,:)];
% b = [1 2 0]';
% Axb = A\b
% Bxb = B\b

%% Problem 1a

a1(1) = 0; % Matrix A Jacobi
a2(1) = 0; % Matrix A Jacobi
a3(1) = 0; % Matrix A Jacobi

b1(1) = 0; % Matrix A G-S
b2(1) = 0; % Matrix A G-S
b3(1) = 0; % Matrix A G-S

c1(1) = 0; % Matrix B Jacobi
c2(1) = 0; % Matrix B Jacobi
c3(1) = 0; % Matrix B Jacobi

d1(1) = 0; % Matrix C G-S
d2(1) = 0; % Matrix C G-S
d3(1) = 0; % Matrix C G-S

% Initial empty matrices for answers
[A1,A2,A3,A4] = deal([]);

for j=2:7
    % Jacobi A
    a1(j) = (1 - 2*a2(j-1) - a3(j-1))/8;
    a2(j) = (2 + a1(j) - a3(j-1))/5;
    a3(j) = (-a1(j) - a2(j))/4;

    A1 = [A1 [a1(j); a2(j); a3(j)]];

    % G-S A
    b1(j) = (1 - 2*b2(j-1) - b3(j-1))/8;
    b2(j) = (2 + b1(j-1) - b3(j-1))/5;
    b3(j) = (-b1(j-1) - b2(j-1))/4;

    A2 = [A2 [b1(j); b2(j); b3(j)]];

    % Jacobi B
    c1(j) = 1 - c2(j-1) - 4*c3(j-1);
    c2(j) = (2 + c1(j) - c3(j-1))/5;
    c3(j) = -8*c1(j) - 2*c2(j);

    A3 = [A3 [c1(j); c2(j); c3(j)]];

    % G-S B
    d1(j) = 1 - d2(j-1) - 4*d3(j-1);
    d2(j) = (2 + d1(j-1) - d3(j-1))/5;
    d3(j) = -8*d1(j-1) - 2*d2(j-1);

    A4 = [A4 [d1(j); d2(j); d3(j)]];

end

% debug displays
% A1
% A2
% A3
% A4

%% Problem 2a

load yalefaces.mat

% X matrix for part C
C = [X(:,37) X(:,532) X(:,1713)];

A5 = C.'*C;

% Debug display
% A5

%% Problem 2b

[V,D] = eigs(A5,3,'LM'); % arrange eigenvalues/vectors largest mag
A6 = [D(1,1) D(2,2) D(3,3)];
A7 = abs(V);

% Debug display
% A6
% A7

%% Problem 2c

[U,S,V] = svd(X);

A8 = abs([U(:,5).'*X(:,37) U(:,5).'*X(:,532) U(:,5).'*X(:,1713)]);

% Debug display
% A8