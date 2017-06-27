function [Data] = load_dataset_interface(foldername)

%%%%% load number of input networks
prompt = 'How many networks do you want to align? ';
num_net = input(prompt); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
networkFolderName = [foldername, '/Net/'];
names={};
for i = 1 : num_net
    prompt = sprintf('Name of %d network:',i);
    Data.Network{i}.name = input(prompt,'s');
    names = [names,Data.Network{i}.name];
    [Data.Network{i}.G,Data.Network{i}.Map] = load_network([networkFolderName, Data.Network{i}.name, '.net']);
    sizes(i) = size(Data.Network{i}.G, 1);
end
Data.Map = names;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the similarity functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
similarityFolderName = [foldername, '/similarity/'];
temp = dir(similarityFolderName);
counter = 1;
for i = 1 : num_net
    for j = i+1 : num_net
        Data.Similarity{counter}.name = [i,'-',j];
        Data.Similarity{counter}.name = [Data.Network{i}.name,'-',Data.Network{j}.name,'.sim'];
        Data.Similarity{counter}.Score = load_similarity([similarityFolderName, Data.Similarity{counter}.name ], sizes(i), sizes(j), Data.Network{i}.Map, Data.Network{j}.Map);
        counter = counter + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G, nodes] = load_network(filename)

f_id = fopen(filename, 'r');
nodes = {};
while 1
    tline = fgetl(f_id);
    if tline == -1
        break;
    end
    temp = strsplit(tline);
    nodes = [nodes,temp{1},temp{2}];
end
fclose(f_id);
nodes = unique(nodes);
%
edges = [];
f_id = fopen(filename, 'r');
while 1
    tline = fgetl(f_id);
    if tline == -1
        break;
    end
    temp = strsplit(tline);
    [~,id1] = ismember(temp{1}, nodes); 
    [~,id2] = ismember(temp{2}, nodes);
    ids=[id1,id2];
    edges = [edges, ids'];
end
fclose(f_id);
numV = max(max(edges));
G = sparse(edges(1,:), edges(2,:), ones(1, size(edges,2)), numV, numV);
G = G + G';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load similarity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Score] = load_similarity(filename, nS, nT, mapS, mapT)
%
edges = [];
f_id = fopen(filename, 'r');
while 1
    tline = fgetl(f_id);
    if tline == -1
        break;
    end
    temp = strsplit(tline,'\t');
    [~,id1] = ismember(temp{1}, mapS); 
    [~,id2] = ismember(temp{2}, mapT);
    if (id1~=0) && (id2~=0)
        ids=[id1,id2,str2num(temp{3})];
        edges = [edges, ids'];
    end
end
fclose(f_id);
Score = sparse(edges(1,:), edges(2,:), edges(3,:), nS, nT);