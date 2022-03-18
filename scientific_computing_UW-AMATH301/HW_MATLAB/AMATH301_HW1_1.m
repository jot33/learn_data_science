% AMATH 301
% HW 1
% University of Washington
% Due 10/7/2021

clear all; close all; clc

%%%%%%%%
% Part 1
%%%%%%%%

% abserr(x) defined at end of script

A1 = abserr(1000)
A2 = abserr(10000)
A3 = abserr(100000)

%%%%%%%%
% Part 2
%%%%%%%%

% Define matrices
A = [1 2; -1 1];
B = eye(2)*2;
C = [[2; 0] [0; 0] [-3; -1]];
D = [1 2; 2 3; -1 0];
x = [1; 0];
y = [0; 1];
z = D(:,1);

A4 = A + B
A5 = 3*x - 4*y
A6 = A*x
A7 = B*(x - y)
A8 = D*x
A9 = D*y + z
A10 = A*B
A11 = B*C
A12 = C*D

%%%%%%%%
% Part 3
%%%%%%%%

p = [0.8 1.5 2.8 3.2 3.5 3.65];
APart3 = ones(1,6)*0.5;

for j=2:50
    APart3 = [APart3; zeros(1,6)];
    for k=1:length(p)
        xn = APart3(size(APart3,1)-1,k);
        APart3(size(APart3,1),k) = p(k)*xn*(1-xn);
    end
end

A13 = APart3(:,1);
A14 = APart3(:,2);
A15 = APart3(:,3);
A16 = APart3(:,4);
A17 = APart3(:,5);
A18 = APart3(:,6);

% Labeled array for debug (comment/uncomment for debug)

% i = reshape([1:50],50,1)    % iteration count
% debugA = [i APart3];
% labels = ["x" "A13(x)" "A14(x)" "A15(x)" "A16(x)" "A17(x)" "A18(x)"];
% debugA = [labels; string(debugA)]


%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%

% Part 1
function y = abserr(x)
    for j=1:x*10
        x = x - 0.1;
    end
    y = abs(0 - x); % output is abs difference between expected and approx.
end