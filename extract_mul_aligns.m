function [Clusters, Map] = extract_mul_aligns(X, offs, eps)
%
dim = offs(length(offs));
networkIds = zeros(1, dim);
for i = 1:(length(offs)-1)
    networkIds((offs(i)+1):offs(i+1)) = i;
end

d = sum(X)';
flags = zeros(1, dim);
Clusters = zeros(1, dim);

[s, order] = sort(-d);

headId = 1;
clusterId = 0;

while 1
    while headId <= dim
        if flags(order(headId)) == 0
            break;
        end
        headId = headId + 1;
    end
    if headId > dim
        break;
    end
    
    rootId = order(headId);
    headId = headId + 1;
    ids = find(X(rootId, :) > eps & flags == 0);
    ids = ids(find(ids ~= rootId));
    if length(ids) > 0
        b = X(ids, rootId);
        A = X(ids, ids);
        A = A - diag(diag(A));
        [x_ids] = mrf_inference(A, b, networkIds(ids));
        scores = b(x_ids) + sum(A(x_ids, x_ids))';
        x_ids = x_ids(find(scores > 0.25));
        if length(x_ids) > 0
            vIds = [rootId, ids(x_ids)];
            clusterId = clusterId + 1;
            Clusters(vIds) = clusterId;
            flags(vIds) = 1;
        end
    else
        h = 10;
    end
end

headId = 1;
while 1
    while headId <= dim
        if flags(order(headId)) == 0
            break;
        end
        headId = headId + 1;
    end
    if headId > dim
        break;
    end
    
    rootId = order(headId);
    headId = headId + 1;
    ids = find(X(rootId, :) > eps & flags == 0);
    ids = ids(find(ids ~= rootId));
    if length(ids) > 0
        b = X(ids, rootId);
        A = X(ids, ids);
        A = A - diag(diag(A));
        [x_ids] = mrf_inference(A, b, networkIds(ids));
        scores = b(x_ids) + sum(A(x_ids, x_ids))';
        x_ids = x_ids(find(scores > eps));
        if length(x_ids) > 0
            vIds = [rootId, ids(x_ids)];
            clusterId = clusterId + 1;
            Clusters(vIds) = clusterId;
            flags(vIds) = 1;
        end
    else
        h = 10;
    end
end

for i = 1:max(Clusters) 
  for j = 1:5
    sub = Clusters((offs(j)+1):offs(j+1));
    off = find(sub == i);
    if length(off) > 0
      Map(i,j) = off(1);
    end
  end
end

function [x_ids] = mrf_inference(A, b, labelIds)
%
nO = 1;
offs = [0];
prevId = labelIds(1);
for off = 2:length(labelIds)
    if labelIds(off) ~= prevId
        nO = nO + 1;
        offs(nO) = off-1;
        prevId = labelIds(off);
    end
end
offs(nO + 1) = length(labelIds);

dim = length(b);
x = ones(dim, 1);
for i = 1:10
    x = normalize_L2(x, offs);
    x = b + A*x;
end
x = normalize_L2(x, offs);
for i = 1:(length(offs)-1)
    [s, tp] = max(x((offs(i)+1):offs(i+1)));
    x_ids(i) = offs(i) + tp;
end

function [x_nor] = normalize_L2(x, offs)
x_nor = x;
for i = 1:(length(offs)-1)
    ids = (offs(i)+1):offs(i+1);
    x_nor(ids) = x_nor(ids)/sum(x_nor(ids));
end