function [X] = match_lift2(W, offs, Para)
% This function implements a pair-wise shape matcher for improving partial
% maps among a collection of objects

% Here we assume the number of underlying points is not fixed

% The input matrix is the scoring matrix between points on all the objects
% offs specify the association between the points and the objects
% Input argument:
%   W: the score matrix between points
%   offs: specify the structure of the block
% Output argument:
%   X: the indicator matrix for the correspondences

dim = offs(length(offs)) + length(offs);
% the linear objective function
C = zeros(dim, dim);
IDX = zeros(1, dim);
for i = 1:(length(offs))
    offs(i) = offs(i) + i;
end
for i = 1:(length(offs)-1)
    IDX((offs(i)+1):(offs(i+1)-1)) = 1;
end
IDX = find(IDX);
C(IDX, IDX) = 0.5 - W;
C = (C+C')/2;
clear W;

% Generate the constraints
[A, b, ids_geq_0] = gen_constraints(offs);


%Variables 
X = zeros(dim, dim);
S = X;

% Perform ADMM to solve the SDP problem
mu = Para.mu_init;
for iter = 1:Para.numIterations
    % Optimize y
    TP = S + mu*X - C;
    y = (A*A')\(A*reshape(TP, [dim*dim,1]) - mu*b);
%    y = TP(ids_A) - mu*b;
    
    % Optimize z_1
    z_1 = -TP(ids_geq_0);
    z_1(find(z_1 < 0)) = 0;
    
    % Optimize S and X
    TP = reshape(A'*y, [dim, dim]);
    TP(ids_geq_0) = TP(ids_geq_0) - z_1;
%    d1 = diag(TP);
    TP = TP + TP' - diag(diag(TP));
%    ids = (0:(dim-1))*dim + (1:dim);
%    TP(ids') = TP(ids') - d1;
    TP = TP + C - mu*X;
    TP = (TP+TP')/2;
    [S, X] = eig(TP);
    x = diag(X)';
    ids = find(x > 0);
    S = S(:,ids)*X(ids,ids)*S(:,ids)';
    X = (S -TP)/mu;
    X = (X+X')/2;
    if mod(iter, 1) == 0
        r = norm(A*reshape(X, [dim*dim,1]) - b);
        fprintf('Finished itertion %d, norm(A*X-b) = %f.\n', iter, r);
    end
    mu = mu*Para.rho;
end
X = X(IDX, IDX);

% Rounding using SVD
% numImages = length(offs)-1;
% maxRank = offs(2:(numImages+1)) - off(1:numImages);
% [U,V] = eig(X);
% eigVals = diag(V)';
% [eigVals, ids] = sort(eigVals);
% ids = ids((length(ids)-maxRank+1):length(ids));
% X = U(:,ids)*V(ids,ids)*U(:,ids)';


function [A, b, ids_geq_0] = gen_constraints(offs)
% Generate the constraints
% A(X) == b
% X(ids_geq_0) >= 0

numNetworks = length(offs)-1;
dim = offs(numNetworks+1);

F = ones(dim, dim);

for i = 1:numNetworks
    ids = (offs(i)+1):offs(i+1);
    F(ids,ids) = 2;
end

ids = 2:dim;
F((ids-1)*dim + 1) = 3;
F((ids-1)*dim + ids) = 3;
F(1,1) = 4;

for i = 1:numNetworks
    F(1, offs(i+1)) = 4;
    F(offs(i+1), offs(i+1)) = 4;
end

for i = 2:dim
    F(i, 1:(i-1)) = 0;
end
F = sparse(F);

% ids_geq_zero
[rows, cols, vals] = find(F == 1);
ids_geq_0 = (cols-1)*dim + rows;

% ids_A
[rows, cols, vals] = find(F == 2);
ids_A = (cols'-1)*dim + rows';
b = zeros(length(cols),1);

[rows, cols, vals] = find(F == 3);
ids_A = [ids_A, (cols'-1)*dim + rows'];
b = [b; ones(length(cols),1)];

% picking up more constraints
% X(1,1) == X(1, offs(i+1))
% X(1,1) == X(offs(i+1), offs(i+1))
rows_A = [];
cols_A = [];
vals_A = [];
b_add = [];
rowId = 0;
for i = 1:numNetworks
    rowId = rowId + 1;
    rows_A = [rows_A, [rowId, rowId]];
    cols_A = [cols_A, [1, (offs(i+1)-1)*dim+1]];
    vals_A = [vals_A, [1,-1]];
    b_add = [b_add; (offs(i+1) - offs(i)-1)];
    rowId = rowId + 1;
    rows_A = [rows_A, [rowId, rowId]];
    cols_A = [cols_A, [1,(offs(i+1)-1)*dim+offs(i+1)]];
    vals_A = [vals_A, [1,-1]];
    b_add = [b_add; (offs(i+1) - offs(i)-1)];
end
A = sparse(1:length(ids_A), ids_A, ones(1,length(ids_A)), length(ids_A), dim*dim);
A_add = sparse(rows_A, cols_A, vals_A, rowId, dim*dim);
A = [A;A_add];
b = [b;b_add];