function [matches] = ConvexAlign(foldername, Para)

%load data
Data = load_dataset_interface(foldername);
if (size(Data.Network(),2)==2)
    pairwise_interface(foldername, Data, Para);
else
    joint_interface(foldername,Data, Para);
end


