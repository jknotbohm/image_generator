% This function obtains the intensity kernel for a selected x-y plane at a given z plane 
% INPUT: z-position, scope's NA, light wavelength, refractive index of objective's immersion
% medium (n)
% OUTPUT: intensity kernel for a selected x-y plane at a given z plane 
% written by Mainak Sarkar, University of Wisconsin-Madison

function [I_xy] = XY_inten_kernel(z, lx, ly, NA, lamda, n)

%% Fix the domain dimensions of XY plane:  (We arrived at these values after a robust trial and error procedure)
  lx = 10 ;
  ly = 10 ;

% generate (x,y) grid on that XY plane:
[x, y] = meshgrid(linspace(-lx/2, lx/2, 100), linspace(-ly/2, ly/2, 100)) ; % # of grid points are achieved by trial and error.

% set v:
v = (2*pi/lamda) * NA * sqrt(x.^2 + y.^2) ;
% set u, which is a constant for a specific XY plane:
u = (2*pi/lamda) * (NA^2/n) * z * ones(size(v));

%% load I_data in a MAT file:
load intensity_kernel.mat I_data

%% Obtain interpolated data set:
F = scatteredInterpolant(I_data(:,1), I_data(:,2), I_data(:,3)) ;

%% Obtain the final intensity kernel for the XY plane under consideration:
I_xy = F(u,v) ;

%% Obtain the image:
I_xy_16bit = uint16(I_xy * 2^16) ; % convert to 16 bit image
