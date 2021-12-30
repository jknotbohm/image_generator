% This function obtains the node coordinates, fibers data, nodal displacements in deformed
% configuration at a given desired step suitable for use in image generator algorithm.
% INPUT: step number, reference nodes and fibers data, displacement data at each steps as extracted from ODB
% OUTPUT: node coordinates, fibers data, nodal displacements in deformed
% configuration at a given desired step suitable for use in image generator algorithm.
% written by Mainak Sarkar, University of Wisconsin-Madison

function [nodes_def, fibers_def, u1, u2, u3] = def_20data(nodes, fibers, get_current, u11, u12, u13, u14, u15, u16, u17, u18, u19, u110, u111, u112, u113, u114, u115, u116, u117, u118, u119, u120, u21, u22, u23, u24, u25, u26, u27, u28, u29, u210, u211, u212, u213, u214, u215, u216, u217, u218, u219, u220, u31, u32, u33, u34, u35, u36, u37, u38, u39, u310, u311, u312, u313, u314, u315, u316, u317, u318, u319, u320)

% [u1, u2, u3, nodes] = inc_def_coords_rev(nodes, u11, u12, u13, u14, u15, u16, u17, u18, u19, u110, u21, u22, u23, u24, u25, u26, u27, u28, u29, u210, u31, u32, u33, u34, u35, u36, u37, u38, u39, u310, get_current) ;
[u1, u2, u3, nodes] = inc_20def_coords_rev(nodes, u11, u12, u13, u14, u15, u16, u17, u18, u19, u110, u111, u112, u113, u114, u115, u116, u117, u118, u119, u120, u21, u22, u23, u24, u25, u26, u27, u28, u29, u210, u211, u212, u213, u214, u215, u216, u217, u218, u219, u220, u31, u32, u33, u34, u35, u36, u37, u38, u39, u310, u311, u312, u313, u314, u315, u316, u317, u318, u319, u320, get_current) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(fibers)
    midpt = .5*(nodes(fibers(j,2),2:4)+nodes(fibers(j,3),2:4));
    midpoint_nodes(j,:) = [length(nodes)+j midpt];
end
% Redo the elements to include the midpoint nodes
fibers = [fibers(:,1) fibers(:,2) midpoint_nodes(:,1)];
g = size(fibers,1)+1 ;
h = 2*size(fibers,1) ;
fibers2 = [(g:h)' midpoint_nodes(:,1) fibers(:,3)];
fibers = [fibers ; fibers2] ;
% Add midpoint nodes to final_nodes
nodes = [nodes; midpoint_nodes];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes_def = nodes ;
fibers_def = fibers ;
