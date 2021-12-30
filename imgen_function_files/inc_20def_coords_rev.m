% This function is a component of def_20data.m
% Written by Mainak Sarkar, UW-Madison

function [u1, u2, u3, nodes_def] = inc_20def_coords_rev(nodes_final, u11, u12, u13, u14, u15, u16, u17, u18, u19, u110, u111, u112, u113, u114, u115, u116, u117, u118, u119, u120, u21, u22, u23, u24, u25, u26, u27, u28, u29, u210, u211, u212, u213, u214, u215, u216, u217, u218, u219, u220, u31, u32, u33, u34, u35, u36, u37, u38, u39, u310, u311, u312, u313, u314, u315, u316, u317, u318, u319, u320, inc) 

%% extract the original undeformed node coordinates:
x = nodes_final(:,2) ;
y = nodes_final(:,3) ; 
z = nodes_final(:,4) ;

U1 = [u11; u12; u13; u14; u15; u16; u17; u18; u19; u110; u111; u112; u113; u114; u115; u116; u117; u118; u119; u120] ;
U2 = [u21; u22; u23; u24; u25; u26; u27; u28; u29; u210; u211; u212; u213; u214; u215; u216; u217; u218; u219; u220] ;
U3 = [u31; u32; u33; u34; u35; u36; u37; u38; u39; u310; u311; u312; u313; u314; u315; u316; u317; u318; u319; u320] ;


for jj = 1:(inc+1)
    if jj == 1
u1 = U1(1,:) ;
u2 = U2(1,:) ;
u3 = U3(1,:) ;
    else
u1 = 0 + U1(jj,:) ;
u2 = 0 + U2(jj,:) ;
u3 = 0 + U3(jj,:) ;    % NO ADDITION
    end
end

%% Report non-uniform deformed node coordinates (fiber end points / nodes)
for i = 1 : size(nodes_final,1)
x_def(i) = x(i) + u1(i) ;
y_def(i) = y(i) + u2(i) ;
z_def(i) = z(i) + u3(i) ;
nodes_def(i,:) = [i x_def(i) y_def(i) z_def(i)] ;
end

