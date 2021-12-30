% This function generates deformed nodes and fibers for the artificial
% rigid body translation. 
% INPUT: nodes, fibers, amount of artificial translation prescribed
% OUTPUT: translated nodes and fibers 
% Written by Mainak Sarkar, University of Wisconsin-Madison

function [nodes_def,fibers_def] = ata_def_nodes(nodes,fibers,x,y) 

for i = 1:size(nodes,1)
    nodes(i,2) = nodes(i,2) + x ;
    nodes(i,3) = nodes(i,3) + y ;
end

nodes_def = nodes ;

fibers_def = fibers ;

