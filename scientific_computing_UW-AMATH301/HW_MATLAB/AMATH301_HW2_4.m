% Setup for A1,A2
clear all; close all; clc

% set given variables
[r1,r2,r3,r4,r5,r6,V1,v2,v3] = deal(20,10,25,10,30,40,0:2:100,0,200);

% setup A matrix based on expanded equations
A = [r1+r2+r6 -r1 -r2; -r1 r1+r3+r4 -r4; -r2 -r4 r2+r4+r5]
[L,U,P] = lu(A);
A1 = zeros(length(V1),3);
A2 = A1;

% solve for currents using backslash and LU methods
for v1=V1
    k = v1/2 + 1;
    b = [v1;v2;v3];
    A1(k,:) = A\b;
    A2(k,:) = U\(L\(P*b)); % x = U\y and y = L\(P*b)
end
A1
A2

% Setup for A3-A5
A3 = zeros(size(A1)); % Initial Jacobi array
A4 = A3; % Initial G-S array
A5 = [0 0]; % Initial iterative success counts

tol = 10e-6; % iterative success tolerance

% initialize currents to zero for iterative method
i1(1)=0;
i2(1)=0;
i3(1)=0;

% loop through using Jacobi iterations
for v1=V1
    k = v1/2 +1;
    i1(k+1) = (v1 + i2(k)*r1 + i3(k)*r2)/(r1 + r2 + r6);
    i2(k+1) = (v2 + i3(k)*r4 + i1(k)*r1)/(r1 + r3 + r4);
    i3(k+1) = (v3 + i2(k)*r4 + i1(k)*r2)/(r2 + r4 + r5);
    A3(k,:) = [i1(k) i2(k) i3(k)];

    A5(1) = k; % update iteration count

    ierr = [i1(k+1)-i1(k) i2(k+1)-i2(k) i3(k+1)-i3(k)];
    if norm(ierr) < tol
        break
    end
end

% loop through using G-s iterations
for v1=V1 % V1=1:2:100, v1 = current iteration of V1
    k = v1/2 + 1; % current iteration

    i1(k+1) = (v1 + i2(k)*r1 + i3(k)*r2)/(r1 + r2 + r6);
    i2(k+1) = (v2 + i3(k)*r4 + i1(k+1)*r1)/(r1 + r3 + r4);
    i3(k+1) = (v3 + i2(k+1)*r4 + i1(k+1)*r2)/(r2 + r4 + r5);

    A4(k,:) = [i1(k) i2(k) i3(k)];
    A5(2) = k; % update iteration count

    ierr = [i1(k+1)-i1(k) i2(k+1)-i2(k) i3(k+1)-i3(k)];

    if norm(ierr) < 10e-6
        break
    end
end

A3
A4
A5

norm(A4(51,:) - A4(50,:))