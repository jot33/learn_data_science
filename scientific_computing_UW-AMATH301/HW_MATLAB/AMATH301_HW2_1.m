% AMATH 301
% HW 2
% University of Washington
% Due 10/14/2021

clear all; close all; clc
%% Problem 1a,b - Backslash and LU solve methods

% Define resistance and voltage values
[r1,r2,r3,r4,r5,r6,v2,v3] = deal(20,10,25,10,30,40,0,200);
V1 = 0:2:100;   % uppercase V1 for array of inputs to v1 variable

% Solve equations symbolically and store coefficients in arrays
syms I1 I2 I3;
eqn1 = r6*I1 + r1*(I1-I2) + r2*(I1-I3);
[C1,~] =coeffs(eqn1);
eqn2 = r3*I2 + r4*(I2-I3) + r1*(I2-I1);
[C2,~] = coeffs(eqn2);
eqn3 = r5*I3 + r4*(I3-I2) + r2*(I3-I1);
[C3,~] = coeffs(eqn3,[I1 I2 I3]);

% Setup matrices and arrays for backslash and LU solve methods
A = [C1;C2;C3] % A matrix for A\b function
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
% x = 1:51;
% plot(x,A1(:,1)); hold on
% plot(x,A1(:,2));
% plot(x,A1(:,3)); hold off

% DEBUG DISPLAY - A2 Output
% disp(double(A2))
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

% Solve equations with respect to variables I1,I2,I3
syms v1; % lower case v1 for single variable
i1 = solve(eqn1 - v1,I1);
i2 = solve(eqn2 - v2,I2);
i3 = solve(eqn3 - v3,I3);

% DEBUG DISPLAY - Symbolic solve
disp([i1;i2;i3])

% Extract coefficients
[i1,var1] = coeffs(i1);
[i2,var2] = coeffs(i2);
[i3,var3] = coeffs(i3);

% DEBUG DISPLAY - Coeff extraction
disp([i1;var1])
disp([i2;var2])
disp([i3;var3])

%% Problem 1b - Jacobian Solver

% Initialize jaocbian arrays
i1j(1) = 0;
i2j(1) = 0;
i3j(1) = 0;

% DEBUG DISPLAY - Iteration equations
% i1j(2) = i1(1)*i2j(1)+ i1(2)*i3j(1) + i1(3)*0
% i2j(2) = i2(1)*i1j(1)+ i2(2)*i3j(1)
% i3j(2) = i3(1)*i2j(1) + i3(2)

%%% Jacobian Solver
for v1 = V1
    k = (v1/2) + 1; % iteration counter
    % Jacobian iterations
    % i1(n) is the nth coefficient of solved i1 eqn
    i1j(k+1) = i1(1)*i2j(k)+ i1(2)*i3j(k) + i1(3)*v1;
    i2j(k+1) = i2(1)*i1j(k)+ i2(2)*i3j(k);
    i3j(k+1) = i3(1)*i2j(k) + i3(2);
    A3(k,:) = [i1j(k) i2j(k) i3j(k)];
    
    if abs(i1j(k+1)-i1j(k)) && ... % check all variables against tol
       abs(i2j(k+1)-i2j(k)) && ...
       abs(i3j(k+1)-i3j(k)) < tol
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