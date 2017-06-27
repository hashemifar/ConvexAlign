function [X_sdp, offs] = joint_multi_convex_align2(Data, Para)
% Jointly aligning multiple graphs
% 
numNetworks = length(Data.Network);
offs = zeros(1, numNetworks+1);
for i = 1 : numNetworks
    offs(i+1) = size(Data.Network{i}.G, 1);
end
%
for i = 1 : numNetworks
    offs(i+1) = offs(i) + offs(i+1);
end
%
for edgeId = 1 : length(Data.Similarity)
    Data.Similarity{edgeId}.Score = normalize_score(Data.Similarity{edgeId}.Score);
end
%
for edgeId = 1:length(Data.Similarity)
    match = Data.Similarity{edgeId};
    [sId, tId] = parse_name(match.name, Data);
    init_matches{edgeId}.sId = sId;
    init_matches{edgeId}.tId = tId;
    init_matches{edgeId}.corres = pairwise_master(...
        Data.Network{sId}.G,...
        Data.Network{tId}.G,...
        match.Score,...
        Para);
end
%%%
% Apply ADMM to solve the sdp formulation
%%%
for outIter = 1 : Para.numOuterIterations
    dim = offs(numNetworks+1);
    X_data = sparse(1:dim, 1:dim, ones(1,dim), dim, dim);
    for eId = 1 : length(init_matches)
        match = init_matches{eId};
        sIds = (offs(match.sId) + 1):offs(match.sId+1);
        tIds = (offs(match.tId) + 1):offs(match.tId+1);
        X_st = sparse(match.corres(1,:),...
            match.corres(2,:),...
            ones(1, size(match.corres,2)),...
            length(sIds),...
            length(tIds));
        X_data(sIds, tIds) = X_st;
        X_data(tIds, sIds) = X_st';
    end
    [eigVecs, eigVals] = eig(full(X_data));
    s = diag(eigVals);
    goodSpecIds = find(s >= numNetworks/2);
    eigVecs2 = eigVecs(:, goodSpecIds);
    eigVals2 = eigVals(goodSpecIds, goodSpecIds);
    X_sdp = eigVecs2*eigVals2*eigVecs2';
    X_sdp = (X_sdp + X_sdp')/2;
        % Now go back to solve the augmented Lagrangian of the pair-wise
    % matching problem
    for edgeId = 1:length(init_matches)
        match = init_matches{edgeId};
        sIds = (offs(match.sId) + 1):offs(match.sId+1);
        tIds = (offs(match.tId) + 1):offs(match.tId+1);
        X_st = X_sdp(sIds, tIds);
        X_st = X_st.*double(X_st > 0.5); % Remove ineffective corres.
        Score = combine_score(Data.Similarity{edgeId}.Score, X_st, Para.lambda_mul);
        corres = pairwise_master(Data.Network{match.sId}.G,...
            Data.Network{match.tId}.G,...
            Score, Para);
        init_matches{edgeId}.corres = corres;
    end
end

function [sId, tId] = parse_name(name, data)
%
temp = strsplit(name,'-');
[~,sId] = ismember(temp{1},data.Map);
temp = strsplit(temp{2},'.');
[~,tId] = ismember(temp{1},data.Map);
%name1 = name(1);
%name2 = name(3);
%sId = int32(name1) - int32('A') + 1;
%tId = int32(name2) - int32('A') + 1;
%
function [Score_out] = normalize_score(Score_in)
%
[ns, nt] = size(Score_in);
[rows, cols, vals] = find(Score_in);
vals = vals/max(vals);
vals = power(max(vals, 0), 0.5);
Score_out = sparse(rows, cols, vals, ns, nt);
Score_out = full(Score_out) + 1e-4;

function [Score] = combine_score(score_in, X_st, alpha)
%
Score = score_in;
[ns,nt] = size(score_in);
Score = reshape(Score, [1, ns*nt]);
X_st = reshape(X_st, [1, ns*nt]);
ids = find(X_st > 0);
Score(ids) = Score(ids) + alpha*X_st(ids);
Score = reshape(Score, [ns, nt]);