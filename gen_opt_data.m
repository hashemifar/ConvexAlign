function [W, offs] = gen_opt_data(Data, IDX)
%
numObjects = length(Data.Network);
offs = zeros(1, numObjects+1);

off = 0;
for i = 1:numObjects
    offs(i+1) = offs(i) + length(IDX{i});
    off = off + length(Data.Network{i}.vertices);
end

W = eye(offs(numObjects+1));

for mapId = 1:length(Data.Maps)
    map = Data.Maps{mapId};
    sId = map.sourceNetworkId;
    tId = map.targetNetworkId;
    X_block = map.G_hub_align(IDX{sId}, IDX{tId});
    sIds = (offs(sId)+1):offs(sId+1);
    tIds = (offs(tId)+1):offs(tId+1);
    W(sIds, tIds) = X_block;
    W(tIds, sIds) = X_block';
end

