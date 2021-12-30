%% This function crops out blank periphery from the image of fiber network (simplified for use on square images)
% INPUT: Uncropped image, crop extent beyond network boundary
% OUTPUT: Cropped image, description of cropping window
% written by Mainak Sarkar, University of Wisconsin-Madison

function [Z_cropped, col1, col2, row1, row2] = crop_network_modified(Z, ex_sp)

d = min(Z(:,1)) ;

for op = 1:size(Z,2)
lg = any(Z(:,op)>d) ;
if lg == 1
break ;
end
end
col1 = op ;

for op = size(Z,2):-1:1
lg = any(Z(:,op)>d) ;
if lg == 1
break ;
end
end
col2 = op ;

for op = 1:size(Z,1)
lg = any(Z(op,:)>d) ;
if lg == 1
break ;
end
end
row1 = op ;

for op = size(Z,1):-1:1
lg = any(Z(op,:)>d) ;
if lg == 1
break ;
end
end
row2 = op ;


% row1 = row1 - ex_sp ;
row1 = ( row2 - (col2 - col1) ) - ex_sp ;
row2 = row2 + ex_sp ;
col1 = col1 - ex_sp ;
col2 = col2 + ex_sp ;

Z_cropped = Z((row2 - (col2 - col1) ) : row2, col1:col2) ; % obtain the cropped image 
