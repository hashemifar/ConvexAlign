The program ConvexAlign finds a global alignment between multiple of networks.   
Please cite: Hashemifar, Somaye, Qixing Huang, and Jinbo Xu. "Joint Alignment of Multiple Protein–Protein Interaction Networks via Convex Optimization." Journal of Computational Biology 23.11 (2016): 903-911.

+++++++++++++++++++++++++++++++++++++++

Format of input files:
•	The network files are identified by the suffix ‘.net’. These files consist of two columns where each row represents an interaction. Column are separated by tab or space.
•	The similarity file of two networks ‘A.net’ and ‘B.net’ is identified by the name ‘A-B.sim’. This file includes three columns where first column, second column and third column respectively represent a protein from first network i.e. A, a protein from the second network i.e. B and their similarity score. Columns are separated by tab or space.
•	Given N network, there should be N(N-1)/2 similarity file for each pair of input networks.
Please notice if the order of input network is A and then B, the name of the corresponding similarity file should be ‘A-B.sim’ not ‘B-A.sim’.

+++++++++++++++++++++++++++++++++++++++

How to run ConvexAlign:
•	Load file parameters.mat
•	Inside convexAlign folder make a work directory with an arbitrary name.  Inside your work directory make two other folders named ‘Net’ and ‘Sim’. Put all the input networks in folder ‘Net’. Put all the pairwise similarities in folder ‘Sim’. 
•	Next run following command:
joint_interface(‘name of folder’, Para);

•	ConvexAlign first asks for the number of input networks that you want to align. Then It asks for name of input networks. While giving the name of networks, do not include the suffix ‘.net’.
•	The output alignment is saved in your work directory under this name: 
                                                   name1-name2- …. -nameK.align 
where name1 to nameK are the names of your K input networks. Output file is a tab separated file in which each line indicates a cluster of matching nodes from input networks. Each cluster at most contains one node from each network. 

+++++++++++++++++++++++++++++++++++++++

There is a sample of three networks in folder data. Below is a run for this sample data:

\>\> joint_interface('data', Para)\;
\nHow many networks do you want to align? 3
\nName of 1 network:A
Name of 2 network:B
Name of 3 network:C

+++++++++++++++++++++++++++++++++++++++

If you have any questions or if you have found any problem in the code, please email Hashemifar@ttic.edu. 

