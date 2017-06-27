function [CORRES] = joint_match_extraction(X, offs, Para)
%
numNetworks = length(offs) - 1;
dim = offs(numNetworks + 1);
for i = 1:dim
    X(i,i) = 0;
end
CORRES = [];
while 1
    corres = extract_match(X, offs, Para.min_val);
    if length(find(corres > 0)) < 2
        break;
    end
    CORRES = [CORRES; corres];
    for i = 1 : length(corres)
        if corres(i) > 0
            corres(i) = offs(i) + corres(i);
        end
    end
    corres = corres(find(corres >0));
    X(corres, :) = 0;
    X(:, corres) = 0;
end

function [corres] = extract_match(X, offs, min_val)
%
corres = -ones(1,length(offs)-1);
[s,rootId] = max(sum(X));
if max(X(rootId,:)) < min_val
    return;
end
ids = [rootId];
i1 = find_offsets(ids(1), offs);
corres(i1) = ids(1) - offs(i1);
while 1
    Y = X(ids, :);
    for j = 1:length(corres)
        if corres(j) > 0
            Y(:, (offs(j)+1):offs(j+1)) = 0;
        end
    end
    if size(Y,1) > 1
        [s,nextId] = max(mean(Y));
    else
        [s,nextId] = max(Y);
    end
    if s < min_val
        break;
    end
    ids = [ids, nextId];
    len = length(ids);
    i2 = find_offsets(ids(len), offs);
    corres(i2) = ids(len) - offs(i2);
end

function [index] = find_offsets(id, offs)
for i = 1:(length(offs)-1)
    if offs(i) < id && offs(i+1) >= id
        index = i;
        break;
    end
end
h = 10;