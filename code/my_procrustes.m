function [d, Z] = my_procrustes(X, Y)
%   [d, Z] = my_procrustes(X, Y) computes the optimal translation, isomorphic
%           scaling and rotation to align Y to X using the following steps:
%     1) Center each shape: Xc = X - mean(X), Yc = Y - mean(Y)
%     2) Scale each shape: Xn = Xc / ||Xc||, Yn = Yc / ||Yc||
%     3) Compute rotation: A = Xn' * Yn, then via SVD: [U,S,V] = svd(A)
%                               and Q = V * U'
%     4) Compute optimal scale:  b = (trace(S) * ||Xc||) / ||Yc||
%     5) Compute aligned shape:  Z = ||Xc||*trace(S) * Yn * Q + mean(X)
%     6) Compute Procrustes distance: d = 1 - (trace(S))^2
%
%   Inputs:
%       X - An N x p matrix (reference shape)
%       Y - An N x p matrix (target shape to be aligned to X)
%
%   Outputs:
%       d - Normalized residual error (Procrustes distance)
%       Z - Transformed version of Y
%

%% 2. Align the shapes to the approximate mean shape
%% a. Calculate the centroid of each shape (or set of landmarks).
muX = mean(X, 1);
muY = mean(Y, 1);

%% b. Align all shapes centroid to the origin.
Xc = X - muX;
Yc = Y - muY;

%% c. Normalize each shapes centroid size
% Centered shape is normalized by its Frobenius norm 
normX = sqrt(sum(sum(Xc.^2)));
normY = sqrt(sum(sum(Yc.^2)));

% Normalize the shapes (isomorphic scaling)
Xn = Xc / normX;
Yn = Yc / normY;

%% d. Rotate each shape to align with the newest approximate mean.
% Compute optimal rotation
A = Xn' * Yn; % (p x p) matrix A = Xn' * Yn
[U, S, V] = svd(A); % Singular value decomposition: A = U * S * V' 
Q = V * U'; % Optimal rotation matrix Q = V * U'

% Compute scale factor
b = trace(S) * normX / normY;

% Compute the transformed Y
Z = b * (Yc * Q) + muX;

% Compute normalized residual error, Procrustes distance
d = 1 - trace(S)^2;
end
