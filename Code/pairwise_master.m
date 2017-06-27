function [corres] = pairwise_master(Graph_s, Graph_t, Similarity, Para)
% Para.lambda_edge : The weight associated with the edge-preservation term
% Para.flag_fast: Whether we try a fast method instead of admm (by default
%                (1: fast, 0: admm)
% 
if isfield(Para, 'flag_fast') == 1 || Para.flag_fast == 1
    X_st = mrf_align_admm2(Graph_s, Graph_t, Similarity,...
        Para.lambda_edge,...
        Para.mu);
else
    Data = mrf_align_opt_data(Graph_s, Graph_t, Similarity);
    X_st = mrf_align_admm(Data, Para.lambda_edge);
end

% Round the fractional solution into an integer solution
ns = size(Graph_s, 1);
nt = size(Graph_t, 1);
rowsol = Hungarian(1-X_st);
corres = [1:ns; rowsol];
