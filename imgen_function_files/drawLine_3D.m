% This function draws a line joining two nodes in 3D space
% INPUT: End nodes
% OUTPUT: 3D line plotted
% written by Mainak Sarkar, University of Wisconsin-Madison

function [] = drawLine_3D(p1, p2)

pts = [p1; p2];
plot3(pts(:,1), pts(:,2), pts(:,3), 'white', 'LineWidth', 4.25e-4)    % one point is 1/72 inch. 150 nm / (1/72) inch ; half-pixel width of fiber almost 150 nm

