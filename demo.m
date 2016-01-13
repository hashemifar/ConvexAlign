function [Clusters, Maps] = demo(Data, Para)
%
[W, offs] = gen_opt_data(Data);
X = match_lift2(W, offs, Para);
[Clusters, Map] = extract_mul_aligns(X, offs, 0.01);