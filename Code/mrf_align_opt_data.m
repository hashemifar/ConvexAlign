function [Data] = mrf_align_opt_data(G_s, G_t, Blast)
% Generate the linear programming data for graph matching
% Input arguments:
%       G_s: the source input graph
%       G_t: the target input graph
%     Blast: the computed pair-wise blast scores

G_s = normalize_graph(G_s);
G_t = normalize_graph(G_t);

%
[rows_S, cols_S, vals_S] = find(G_s);
edges_S = [cols_S'; rows_S'];
%
[rows_T, cols_T, vals_T] = find(G_t);
edges_T = [cols_T'; rows_T'];
%
ns = size(G_s, 1);
nt = size(G_t, 1);
[rows, cols, vals] = find(ones(ns, nt));
mask_data = [rows';cols'];
[E1_rowIds, E2_rowIds, rowCorIds1] = BinaryCons(edges_S, edges_T, mask_data, [ns,nt]);
mask_data = [cols';rows'];
[E3_rowIds, E4_rowIds, rowCorIds2] = BinaryCons(edges_T, edges_S, mask_data, [nt,ns]);

dimCor = length(rowCorIds1);
dimY = size(E1_rowIds, 2);
E1 = sparse(E1_rowIds, 1:dimY, ones(1,dimY), dimCor, dimY);
idsF1 = find(sum(E1')>0);
E1 = E1(idsF1, :);
idsF1 = rowCorIds1(idsF1);

E2 = sparse(E2_rowIds, 1:dimY, ones(1,dimY), dimCor, dimY);
idsF2 = find(sum(E2')>0);
E2 = E2(idsF2, :);
idsF2 = rowCorIds1(idsF2);

dimCor = length(rowCorIds2);
dimY = size(E3_rowIds, 2);
E3 = sparse(E3_rowIds, 1:dimY, ones(1,dimY), dimCor, dimY);
idsF3 = find(sum(E3')>0);
E3 = E3(idsF3, :);
idsF3 = rowCorIds2(idsF3);

E4 = sparse(E4_rowIds, 1:dimY, ones(1,dimY), dimCor, dimY);
idsF4 = find(sum(E4')>0);
E4 = E4(idsF4, :);
idsF4 = rowCorIds2(idsF4);

TP = reshape(1:(ns*nt), [ns, nt])';
TP = reshape(TP, [1, ns*nt]);
idsF3 = TP(idsF3);
idsF4 = TP(idsF4);

Mask = reshape(1:(ns*nt), [ns, nt]);

[v1Ids, s, v] = find(E1);
[v2Ids, s, v] = find(E2);
[v3Ids, s, v] = find(E3);
[v4Ids, s, v] = find(E4);
v1Ids = Mask(idsF1(v1Ids'));
v2Ids = Mask(idsF2(v2Ids'));
v3Ids = Mask(idsF3(v3Ids'));
v4Ids = Mask(idsF4(v4Ids'));
num = length(v1Ids);

M1 = sparse(min(v1Ids, v2Ids), max(v1Ids, v2Ids), 1:num);
M2 = sparse(min(v3Ids, v4Ids), max(v3Ids, v4Ids), 1:num);
[v1Ids, v2Ids, order12] = find(M1);
[v3Ids, v4Ids, order34] = find(M2);
E1 = E1(:, order12');
E2 = E2(:, order12');
E3 = E3(:, order34');
E4 = E4(:, order34');


Data.E = [E1;E2;E3;E4];
Data.offs = [0, length(idsF1), length(idsF2), length(idsF3), length(idsF4)];
Data.idsF = [idsF1, idsF2, idsF3, idsF4];
for i = 1:4
    Data.offs(i+1) = Data.offs(i) + Data.offs(i+1);
end
Data.w_y = ones(size(Data.E, 2), 1);
Data.W_X = Blast;

function [G_nor] = normalize_graph(G)
% We will high-weights to edges 
d = max(full(sum(G)), 1);
ids_non_hub = find(d < max(d)/12);
[rows, cols, vals]= find(G);
d = d/max(d);
vals_nor = max(d(rows'), d(cols'));
if 1
    vals_nor = sqrt(d(rows').*d(cols'));
end

n = size(G, 1);
G_nor = sparse(rows, cols, vals_nor, n, n);
G_nor(ids_non_hub, ids_non_hub) = 0;
G_nor = sparse(G_nor);