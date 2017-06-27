function make_output(Data, matches, foldername)

name=[Data.Network{1}.name];
for i = 2:size(Data.Network,2)
    name=[name,'-',Data.Network{i}.name]; 
end
name = [foldername,'/',name,'.align'];
%
fileID = fopen(name,'w');
for i = 1:size(matches,1)
    for j = 1:size(Data.Network,2)
        if (matches(i,j)==-1)
            fprintf(fileID,'-\t');
        else
            fprintf(fileID,'%s\t',Data.Network{j}.Map{matches(i,j)});
        end
    end
    fprintf(fileID,'\n');
end