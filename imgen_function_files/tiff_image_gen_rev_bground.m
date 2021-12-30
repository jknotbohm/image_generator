% This function generates all the synthetic images
% INPUT: Input of nodes, fibers of 3D network, wavelength, scope
% parameters, objective details (NA, refractive index of immersion medium),
% z-resulution limits.
% OUTPUT: Uncropped images at each z plane, all x-y kernels at respective z
% planes.
% written by Mainak Sarkar, University of Wisconsin-Madison

function [J, kernel] = tiff_image_gen_rev_bground(get_current, z_const, Z_LIMIT, resolution, nodes, fibers, lx, ly, lz, magnification, NA, FF, GG, lamda, n)  % generates all the tiff images
set(0,'DefaultFigureVisible','on')  % allows figures to display
h = figure(get_current+1) ;

%% develop the 3D model:
counter = 0 ;
counter1 = 0 ;
for i = 1:size(fibers,1)
   for j = 1:size(nodes,1)
   if fibers(i,2)==nodes(j,1)
       counter = counter + 1 ;
       X(counter,:) = [nodes(j,2) nodes(j,3) nodes(j,4)] ;
   end
   end
   for j = 1:size(nodes,1)
   if fibers(i,3)==nodes(j,1) 
       counter1 = counter1 + 1 ;
       Y(counter1,:) = [nodes(j,2) nodes(j,3) nodes(j,4)] ;
   end
   end
drawLine_3D(X(counter,:), Y(counter1,:))   % line fibers
hold on 
end 

%% generate images:
for ij = 1 : size(Z_LIMIT,1)
disp(['start of generation of image sub-slice ',num2str(ij)])
P = Z_LIMIT(ij,1) + 0.000 ; % for avoiding the pixel saturation, this modification (+ 0.000) is enforced
Q = Z_LIMIT(ij,2) ;
xlim([-lx/2 lx/2])
ylim([-ly/2 ly/2])
zlim([P Q]) % slice the 3D domain for 'z' belongs to [P, Q].
view(2) % top view
daspect([magnification magnification 1])  % control the relative lengths of one data unit along each axis (use it / recommended) ** MANDATORY
axis off
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);  % recommended

% minimize figure window 
set(gcf,'WindowState','minimized');

%% IMAGE GENERATION:

set(gcf,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color

% A figure can be converted into an image using the GETFRAME and FRAME2IM functions in MATLAB. 
% For example, the following code creates a snapshot of the current figure and writes it to an image in variable "X1" with colormap stored in "Map":
F = getframe(gcf) ;

[X1, Map] = frame2im(F) ;
imwrite(X1, 'pout_png.png', 'XResolution', resolution, 'YResolution', resolution)
org_img = imread('pout_png.png') ; 

tiff_image_8bit = rgb2gray(org_img); 

J2 = imadjust(tiff_image_8bit, [0,1]) ;

disp('start of convolution operation with respective XY-kernel')
[J{ij}, kernel{ij}] = conv_with_intensity(J2, P, Q, z_const, lx, ly, NA, FF, GG, lamda, n) ;
disp('end of convolution operation with respective XY-kernel')

disp(['end of generation of image sub-slice ',num2str(ij)])

end
