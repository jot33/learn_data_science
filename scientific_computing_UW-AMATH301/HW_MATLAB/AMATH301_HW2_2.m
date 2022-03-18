% AMATH 301
% HW 2
% University of Washington
% Due 10/14/2021

clear all; close all; clc
%% Problem 1a,b - Backslash and LU solve methods

% Define resistance and voltage values
[r1,r2,r3,r4,r5,r6,v2,v3] = deal(20,10,25,10,30,40,0,200);
V1 = 0:2:100;   % uppercase V1 for array of inputs to v1 variable

% Setup matrices and arrays for backslash and LU solve methods
A = [r1+r2+r6 -r1 -r2
    -r1 r1+r3+r4 -r4
    -r2 -r4 r2+r4+r5] % A matrix for A\b function
A1 = zeros(length(V1),3);
A2 = zeros(size(A1)); % initialize 51,3 zero array
[L,U,P] = lu(A); % solve for L,U,P matrices based on previously found A

for v1 = V1 % lowercase v1 is current value from V1 array
    k = v1/2 + 1; % iteration number
    b = [v1;v2;v3];
    
    % Backslash solver
    A1(k,:) = A\b;
    
    % LU Solver
    % y = L\(P*b) and x = U\y => x = U\(L\(P*b))
    A2(k,:) = U\(L\(P*b));
end

% DEBUG DISPLAY - A1 output (symbolic solve gen rational -> cast double)
% disp(double(A1))
% disp(A1)
% x = 1:51;
% plot(x,A1(:,1)); hold on
% plot(x,A1(:,2));
% plot(x,A1(:,3)); hold off

% DEBUG DISPLAY - A2 Output
% disp(double(A2))
% disp(A2)
% x = 1:51;
% plot(x,A2(:,1)); hold on
% plot(x,A2(:,2));
% plot(x,A2(:,3)); hold off

%% Problem 1b,c - Iterative Setup

%%% Iterative Solve setup
tol = 10e-6; % tolerance for iterative solve methods
A3 = zeros(size(A1)); % Jacobian solve initial array
A4 = A3; % G-S solve initial array
A5 = [0 0]; % Initialized array for iteration counts
%% Problem 1b - Jacobian Solver

% Initialize jaocbian arrays
i1j(1) = 0;
i2j(1) = 0;
i3j(1) = 0;

%%% Jacobian Solver
for v1 = V1
    k = v1/2 + 1; % iteration counter
    % Jacobian iterations
    i1j(k+1) = (v1 + r1*i2j(k) + r2*i3j(k))/(r1+r2+r6);
    i2j(k+1) = (v2 + r1*i1j(k) + r4*i3j(k))/(r1+r3+r4);
    i3j(k+1) = (v3 + r2*i1j(k) + r4*i2j(k))/(r2+r4+r5);
    A3(k,:) = [i1j(k) i2j(k) i3j(k)];
    
    i1err = abs(i1j(k+1)-i1j(k));
    i2err = abs(i2j(k+1)-i2j(k));
    i3err = abs(i3j(k+1)-i3j(k));
    if and(and(i1err<tol,i2err<tol),i3err<tol)
        
        A5(1,1) = k; % save iteration count in A5
        break
    end
end

% DEBUG DISPLAY - A4 output
size(A4)
disp(A4)
disp(A5)
%% Problem 1b - Gauss-Seidel Solver

% Initialize Gauss-Seidel (G-S) arrays
i1g(1) = 0;
i2g(1) = 0;
i3g(1) = 0;

%%% Jacobian Solver
for v1 = V1
    k = (v1/2) + 1; % iteration counter
    % Gauss-Seidel iterations
    i1g(k+1) = i1(1)*i2g(k)+ i1(2)*i3g(k) + i1(3)*v1;
    i2g(k+1) = i2(1)*i1g(k+1)+ i2(2)*i3g(k);
    i3g(k+1) = i3(1)*i2g(k+1) + i3(2);
    A3(k,:) = [i1g(k+1) i2g(k+1) i3g(k+1)];
    
    if abs(i1j(k+1)-i1j(k)) && ... % check all variables against tol
       abs(i2j(k+1)-i2j(k)) && ...
       abs(i3j(k+1)-i3j(k)) < tol
        A5(1,2) = k; % save iteration count in A5
        break
    end
end