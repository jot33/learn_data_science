% AMATH 301
% HW 3
% University of Washington
% Due 10/25/2021

clear all; close all; clc

% load yalefaces
load yalefaces.mat

%% 1a compute 100 x 100 correlation matrix

% Extract first 100 images (columns) and find correlation by C = X.'*X
A1 = X(:,1:100).'*X(:,1:100); 
sizeA1 = size(A1) % display size of A1

%% 1b determine highest and lowest correlations

mxcor = 0; % max correlation
for j = 1:length(A1)-1
    for k = j+1:length(A1)
        if A1(j,k) > mxcor
            mxcor = A1(j,k);
            A2 = [j k];
        end
    end
end
mxcor % display max correlation
A2 % display A2

mncor = mxcor; % minimim correlation
for j = 1:length(A1)-1
    for k = j+1:length(A1)
        if A1(j,k) < mncor
            mncor = A1(j,k);
            A3 = [j k];
        end
    end
end

mncor
A3

%% 1c repeat 1a with 10 x 10 correlation matrix

% X matrix for part C
Xc = [X(:,1) X(:,313) X(:,512) X(:,5) ...
    X(:,2400) X(:,113) X(:,1024) X(:,87) X(:,314) X(:,2005)];

A4 = Xc.'*Xc
sizeA4 = size(A4)
%% 1d create Y = XX^T, find 6 largest eigs
Y = X*X.';
sizeY = size(Y)
[V,D] = eigs(Y,6,'LM'); % arrange eigenvalues/vectors largest mag
A5 = abs(V); % absolute value of eigenvectors

%% 1e SVD matrix X, find 6 principal directions

[U,S,V] = svd(X); % singular value decomp
A6 = U(:,1:6); % first 6 col of unitary matrix U are 6 principal dir

%% 1f norm between eig vector v1 and SVD mode u1

A7 = norm(A5(:,1) - abs(A6(:,1)))

%% 1g compute percentage of variance capture by each of first 6 SVD modes
totalvar = sum(S,'all');
A8 = [];
for j = 1:6
    A8 = [A8 100*S(j,j)/totalvar];
end
A8