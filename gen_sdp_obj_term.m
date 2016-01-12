function [W] = gen_sdp_obj_term(Data)
%
D_00 = power(Data.G_00_aff.*Data.D_00_aff, 0.5);
D_11 = power(Data.G_11_aff.*Data.D_11_aff, 0.5);
D_22 = power(Data.G_22_aff.*Data.D_22_aff, 0.5);
D_33 = power(Data.G_33_aff.*Data.D_33_aff, 0.5);
D_44 = power(Data.G_44_aff.*Data.D_44_aff, 0.5);
%
offs = Data.offs;
dim = offs(length(offs));
W = zeros(dim, dim);
%
ids0 = (offs(1)+1):offs(2);
ids1 = (offs(2)+1):offs(3);
ids2 = (offs(3)+1):offs(4);
ids3 = (offs(4)+1):offs(5);
ids4 = (offs(5)+1):offs(6);
%
W_st = expand_weight_mat(Data.G_01, D_00, D_11);
W(ids0, ids1) = W_st;
W(ids1, ids0) = W_st';

W_st = expand_weight_mat(Data.G_02, D_00, D_22);
W(ids0, ids2) = W_st;
W(ids2, ids0) = W_st';

W_st = expand_weight_mat(Data.G_03, D_00, D_33);
W(ids0, ids3) = W_st;
W(ids3, ids0) = W_st';

W_st = expand_weight_mat(Data.G_04, D_00, D_44);
W(ids0, ids4) = W_st;
W(ids4, ids0) = W_st';

W_st = expand_weight_mat(Data.G_12, D_11, D_22);
W(ids1, ids2) = W_st;
W(ids2, ids1) = W_st';

W_st = expand_weight_mat(Data.G_13, D_11, D_33);
W(ids1, ids3) = W_st;
W(ids3, ids1) = W_st';

W_st = expand_weight_mat(Data.G_14, D_11, D_44);
W(ids1, ids4) = W_st;
W(ids4, ids1) = W_st';

W_st = expand_weight_mat(Data.G_23, D_22, D_33);
W(ids2, ids3) = W_st;
W(ids3, ids2) = W_st';

W_st = expand_weight_mat(Data.G_24, D_22, D_44);
W(ids2, ids4) = W_st;
W(ids4, ids2) = W_st';

W_st = expand_weight_mat(Data.G_34, D_33, D_44);
W(ids3, ids4) = W_st;
W(ids4, ids3) = W_st';


function [G_st] = expand_weight_mat(G_st_in, D_ss, D_tt)
%
G_st = (G_st_in*D_tt + D_ss*G_st_in)/2;
