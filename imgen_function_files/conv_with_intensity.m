% This function does the convolution operation with intensity kernel for
% respective x-y image on and across focal plane
% INPUT: Ground truth image data, limits of z for differential 3D volume, z-position of focus, 
% NA, light wavelength, objective refractive index, 
% OUTPUT: Convolved image, rexpective x-y kernel
% written by Mainak Sarkar, University of Wisconsin-Madison

function [out_img, kernel] = conv_with_intensity(in_img, P, Q, z_const, lx, ly, NA, FF, GG, lamda, n) 

z = 0.5 * (P + Q) - z_const ; % focal plane is z = z_const

%% obtain the kernel:

% Develop the respective x-y intensity kernel:
[I_xy] = XY_inten_kernel(z, lx, ly, NA, lamda, n) ;

%% 2D convolution:
kernel = I_xy ;
out_img = conv2(in_img, kernel, 'same') ;
out_img = im2double(out_img) ;

  


