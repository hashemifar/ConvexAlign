function [corres] = pairwise_interface(foldername, Data, Para)
% Para.lambda_edge : The weight associated with the edge-preservation term
% Para.flag_fast: Whether we try a fast method instead of admm (by default
%                (1: fast, 0: admm)
% 
Similarity = normalize_score(Data.Similarity{1}.Score);
corres = pairwise_master(Data.Network{1}.G, Data.Network{2}.G, Similarity, Para);
make_output(Data, corres',foldername)

function [Score_out] = normalize_score(Score_in)
%
[ns, nt] = size(Score_in);
[rows, cols, vals] = find(Score_in);
vals = vals/max(vals);
vals = power(max(vals, 0), 0.5);
Score_out = sparse(rows, cols, vals, ns, nt);
Score_out = full(Score_out) + 1e-4;
