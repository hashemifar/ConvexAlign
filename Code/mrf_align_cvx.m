function [X] = mrf_align_cvx(Data, lambda)
% This function is the accompanying function of mrf_align_admm. It is used
% to check the correctness of mrf_align_admm
%
% This function solves the following linear program
% maximize     <W_X, X > + <w_y, y>
% subject to   y     >= 0
%              E*y <= F*X
%              X* 1  <= 1
%              X'*1  <= 1
%              X     >= 0

[numV_s, numV_t] = size(Data.W_X);
dim_y = length(Data.w_y);

cvx_begin
variable X(numV_s, numV_t);
variable y(dim_y);
maximize ((lambda*Data.w_y)'*y + sum(sum(Data.W_X.*X)))
subject to
X >= 0;
y >= 0;
Data.E*y <= X(Data.idsF');
X*ones(numV_t, 1) <= 1;
X'*ones(numV_s, 1) <= 1;
cvx_end

fprintf(' energy = [%f, %f]\n', Data.w_y'*y, sum(sum(X.*Data.W_X)));