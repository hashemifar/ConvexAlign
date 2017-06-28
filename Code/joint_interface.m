function [matches] = joint_interface(foldername,Data, Para)

% Compute the data-matrix 
[X, offs] = joint_multi_convex_align(Data, Para);
% Round the data matrix into consistent matches
matches = joint_match_extraction(X, offs, Para);
%make alignment
make_output(Data, matches,foldername)

