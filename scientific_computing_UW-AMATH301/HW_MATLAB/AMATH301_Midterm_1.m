% AMATH 301
% Midterm 1
% University of Washington
% Due 10/30/2021

clear all; close all; clc

%% Problem 1

% Define given matrices
A = [8 2 1; -1 5 1;  1 1 4];
B = [A(3,:); A(2,:); A(1,:)];
b = [1 2 0]';

% tests for debugging
% Axb = A\b
% Bxb = B\b
%% Problem 1a

i1AJ(1) = 0; % Matrix A Jacobi
i2AJ(1) = 0; % Matrix A Jacobi
i3AJ(1) = 0; % Matrix A Jacobi

i1AG(1) = 0; % Matrix A G-S
i2AG(1) = 0; % Matrix A G-S
i3AG(1) = 0; % Matrix A G-S

i1BJ(1) = 0; % Matrix B Jacobi
i2BJ(1) = 0; % Matrix B Jacobi
i3BJ(1) = 0; % Matrix B Jacobi

i1BG(1) = 0; % Matrix C G-S
i2BG(1) = 0; % Matrix C G-S
i3BG(1) = 0; % Matrix C G-S

% Initial empty matrices for answers
[A1,A2,A3,A4] = deal([]);

for j=1:7
    % Jacobi A
    i1AJ(j+1) = 1/8 - i2AJ(j)/4 - i3AJ(j)/8;
    i2AJ(j+1) = 2/5 + i1AJ(j+1)/5 - i3AJ(j)/5;
    i3AJ(j+1) = -i1AJ(j+1)/4 - i2AJ(j+1)/4;

    A1 = [A1 [i1AJ(j+1); i2AJ(j+1); i3AJ(j+1)]];

    % G-S A
    i1AG(j+1) = 1/8 - i2AG(j)/4 - i3AG(j)/8;
    i2AG(j+1) = 2/5 + i1AG(j)/5 - i3AG(j)/5;
    i3AG(j+1) = -i1AG(j)/4 - i2AG(j)/4;

    A2 = [A2 [i1AG(j+1); i2AG(j+1); i3AG(j+1)]];

    % Jacobi B
    i1BJ(j+1) = 1 - i2BJ(j) - 4*i3BJ(j);
    i2BJ(j+1) = 2/5 + i1BJ(j+1)/5 - i3BJ(j)/5;
    i3BJ(j+1) = -8*i1BJ(j+1) - 2*i2BJ(j+1);

    A3 = [A3 [i1BJ(j+1); i2BJ(j+1); i3BJ(j+1)]];

    % G-S B
    i1BG(j+1) = 1 - i2BG(j) - 4*i3BG(j);
    i2BG(j+1) = 2/5 + i1BG(j)/5 - i3BG(j)/5;
    i3BG(j+1) = -8*i1BG(j) - 2*i2BG(j);

    A4 = [A4 [i1BG(j+1); i2BG(j+1); i3BG(j+1)]];
end

A1 = A1(:,1:6);
A2 = A2(:,1:6);
A3 = A3(:,1:6);
A4 = A4(:,1:6);

% debug displays
% A1
% A2
% A3
% A4

%% Problem 2a

load yalefaces.mat

% X matrix for part C
A5 = [X(:,37) X(:,532) X(:,1713)];

A5 = A5.'*A5;

% Debug display
% A5

%% Problem 2b

[A7,A6] = eigs(A5,3,'LM'); % arrange eigenvalues/vectors largest mag
A6 = [A6(1,1); A6(2,2); A6(3,3)];

% Debug display
% A6
% A7

%% Problem 2c

[U,S,V] = svd(X);

A8 = [U(:,5).'*X(:,37); U(:,5).'*X(:,532); U(:,5).'*X(:,1713)];

% Debug display
% A8