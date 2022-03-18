% Setup for A1-A5
clear all; close all; clc

% set variables
[r1,r2,r3,r4,r5,r6,V1,v2,v3] = deal(20,10,25,10,30,40,0:2:100,0,200);
tol = 1e-6;                % iterative success tolerance

% set empty array to fill out
A1 = zeros(length(V1),3);   % Backslash method array
A2 = A1;                    % LU method array
A3 = zeros(size(A1));       % Jacobi array
A4 = A3;                    % G-S array

% setup A matrix based on expanded equations
A = [r1+r2+r6 -r1 -r2; -r1 r1+r3+r4 -r4; -r2 -r4 r2+r4+r5];
[L,U,P] = lu(A);

% initialize "guess" for currents
i1(1)=0;
i2(1)=0;
i3(1)=0;

% initial jacobi and G-S iteration counts
jc(1) = 0; % Jacobi count
gc(1) = 0; % G-S count

% Solver
for v1=V1
    j = v1/2 + 1;       % Current voltage iteration count
    b = [v1;v2;v3];     % Voltage vector

    % Solve backslash method
    A1(j,:) = A\b;

    % Solve LU method
    A2(j,:) = U\(L\(P*b)); % x = U\y and y = L\(P*b)

    % Solve Jacobi iterative method
    for k = 1:200
        % Process equations
        i1(k+1) = (v1 + i2(k)*r1 + i3(k)*r2)/(r1 + r2 + r6);
        i2(k+1) = (v2 + i3(k)*r4 + i1(k)*r1)/(r1 + r3 + r4);
        i3(k+1) = (v3 + i2(k)*r4 + i1(k)*r2)/(r2 + r4 + r5);
        
        jc(j) = k; % update iteration count
        ierr = [i1(k+1)-i1(k) i2(k+1)-i2(k) i3(k+1)-i3(k)]; % error array

        if norm(ierr) < tol
            A3(j,:) = [i1(k+1) i2(k+1) i3(k+1)];
            break
        end
    end

    % Solve G-S iterative method
    for k = 1:200
        % Process equations
        i1(k+1) = (v1 + i2(k)*r1 + i3(k)*r2)/(r1 + r2 + r6);
        i2(k+1) = (v2 + i3(k)*r4 + i1(k+1)*r1)/(r1 + r3 + r4);
        i3(k+1) = (v3 + i2(k+1)*r4 + i1(k+1)*r2)/(r2 + r4 + r5);
        
        gc(j) = k; % update iteration count
        ierr = [i1(k+1)-i1(k) i2(k+1)-i2(k) i3(k+1)-i3(k)]; % error array

        if norm(ierr) < tol
            A4(j,:) = [i1(k+1) i2(k+1) i3(k+1)];
            break
        end
    end
end

A5 = [mean(jc) mean(gc)];

% Debug Displays (comment/uncomment below)

% jc(41:51)
% gc(41:51)
% A
% A1
% A2
% A3
% A4
% A5