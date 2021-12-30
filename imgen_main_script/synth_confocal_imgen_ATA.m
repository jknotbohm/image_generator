% microscope like step by step image generator at different NA
% Can form images in any imposed artificial rigid body translation 
% written by Mainak Sarkar, University of Wisconsin-Madison

close all
clear all  
clc; 
tic

%% INPUT PANEL (Input the properties of objective lens and light)
magnification = 20 ; % used in downsampling step at the end / for objective lens
pix_size = 6.5e-0 ; % in microns unit / for camera / used in downsampling step at the end
lamda = 525 * 10^-3 ; % wavelength in micron unit
n = 1 ; % refractive index of the air or oil or water b/w specimen and the objective
z_const = 0 ; % about z = 0 if z_const = 0 / specify the focal plane here (if different).

% This script is designed to work on one single default value of a range of values of NA. 
% for NA = 0.75:-0.05:0.65   % Numerical Aperture of the objective
for NA = 0.7
% specify translation in the final image along x and y: 
% please use distance unit, like microns here:
x = 3 ; % applied x-displacement in um
y = 2.5 ; % applied y-displacement in um
xyk = 10 ;
xx = x ; yy = y ;

%---------------------------------------------------------------

%% Input switch (reference / deformed):
for get_current = [0:1]  % reference and translated images
% for get_current = [0]  % reference
%% load: 
index_inp = 1003 ;
load(['C:\Users\msarkar3\Downloads\DIC_pjt_code_aesthetic_edits_3Dec21\image_generator\imgen_function_files\',...
    '3D_d0p025f20_Network_350x350x5_Seed_',num2str(index_inp),'\'...
    '3D_d0p025f20_Network_350x350x5_Seed_cropped_',num2str(index_inp),'.mat'],'final_nodes','final_fibers', 'lx', 'ly', 'lz') 

nodes = final_nodes ;

fibers = [final_fibers(:,1), final_fibers(:,2), final_fibers(:,4)] ;

ref_nodes = nodes ;

if get_current == 0 
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
end

resolution = 150 ;

% NA = 0.7 ; % Numerical Aperture of the objective / specified at the beginning

if get_current ~= 0
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
[nodes_def,fibers_def] = ata_def_nodes(nodes,fibers,x,y) ;
nodes = nodes_def ;
fibers = fibers_def ;
end

u  = (-10:0.01:10) ;
I = (sin(u/4)./(u/4)).^4 ;  

hh = 0 ;
for int = .9:-.3:0 
a = find(I>int) ; 
st = u(a(1)) ; % lower limit 
last = u(a(end)) ; % upper limit 
const_fact = ((2*pi)/lamda)*(NA^2/n) ;
hh = hh + 1 ;
z_lb(hh) = (1/const_fact) * st  + z_const ;
z_ub(hh) = (1/const_fact) * last + z_const ;
end

% image @ focal plane z = z_const
P(1) = z_lb(1) ;
Q(1) = z_ub(1) ;
i = 1 ;
for g = 1 : (size(z_ub,2)-1)
i = i + 1 ;
P(i) = z_ub(g) ;
Q(i) = z_ub(g+1) ;
end
for g = 1 : (size(z_lb,2)-1)
i = i + 1 ;
P(i) = z_lb(g+1) ;
Q(i) = z_lb(g) ;
end

Z_LIMIT = [P' Q'] ;
bb = (size(Z_LIMIT,1) - 1) / 2 ;
uu = 0 ;
for jj = size(Z_LIMIT,1)-(bb-1) :1: size(Z_LIMIT,1)
uu = uu + 1 ;
PQ(uu,:) = Z_LIMIT(jj,:) ;
end
PQ = flip(PQ) ;
st = size(Z_LIMIT,1)-(bb-1) - 1 ;
Z_LIMIT = [PQ; Z_LIMIT(1:st, :)] ;

FF = z_ub(end) ;
GG = z_lb(end) ;

% Get the uncropped images at different z planes:
[J, kernel] = tiff_image_gen_rev_bground(get_current, z_const, Z_LIMIT, resolution, nodes, fibers, lx, ly, lz, magnification, NA, FF, GG, lamda, n) ; % generates all the PNG images

for k = 1 : size(J,2)
J{k} = rescale(J{k}, 0, (max(max(J{k})) - min(min(J{k})))) ;
end
%% which method to use?
method = 1 ; % specify 1 or 2
%% path / method 1:
% alternate: use of addition:
if method == 1
for h = 1 : (size(J,2)-1)
if h == 1
a = J{h} ;
end
a = (a + J{h+1}) ;
end
Z = a ;
Z = rescale(Z, 0, 1) ;
end

%% method / path 2:
if method == 2 
J = rescaling_inten(J, 0, 0, 0) ;
for h = 1 : (size(J,2)-1)
if h == 1
a = J{h} ;
end
a = (a + J{h+1}) ;
end
Z = a ;
end

% 'crop out' operation
if get_current == 0
ex_sp = 0 ;
% [Z_cropped, col1, col2, row1, row2] = crop_network(Z, ex_sp) ; % You can
% use this if you let figure window to display in full screen in the process
[Z_cropped, col1, col2, row1, row2] = crop_network_modified(Z, ex_sp) ; % recommended / please use this if you are suppressing figures from display
else
Z_cropped = Z(row1:row2, col1:col2) ; % same window for translated or deformed states
end

%% RE-sampling the image to fit to microscope image requirements:
myDPI = get(groot,'ScreenPixelsPerInch') ; % get your screen's pixels per inch
pix_screen = (1/myDPI)*25.4*1000 ; % in microns
d1 = lx*magnification/(pix_size) * (1) ; % convert to computer screen's pixel size from the camera pixel size  
d2 = ly*magnification/(pix_size) * (1) ; % convert to computer screen's pixel size from the camera pixel size 
d1 = ceil( d1 ) ;
d2 = ceil( d2 ) ;
Z_mssmpl = imresize(Z_cropped, [d1, d2]) ;   % if lx, ly in micron unit

%% real image (sample for training):
ref1 = imread('sample_confocal_image.tif', 4) ;
median_inten_ref = double(median(ref1, 'all')) ;
mean_inten_ref = double(mean(ref1, 'all')) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% convert to 16 bit image:
Z_mssmpl1 = uint16(Z_mssmpl * double(max(ref1(:)))) ; % convert to 16 bit image

%% ADD NOISE NOW:
median_inten_curr = median(Z_mssmpl1, 'all') ;
mean_inten_curr = mean(Z_mssmpl1, 'all') ;
mean_noise = (22.285/median_inten_ref) * median_inten_curr ;   % (standard Gaussian, std: 22.285, 10.178 sourced from representative real images)
sd_noise = (10.178/median_inten_ref) * median_inten_curr ;     % 
var_noise = (sd_noise)^2 ;
noised_base = uint16(sqrt(double(var_noise)).*randn(size(Z_mssmpl1)) + double(mean_noise)) ;
Z_mssmpl2 = imadd(Z_mssmpl1, noised_base, 'uint16') ;

%% store final tiff image in file
imwrite(Z_mssmpl2,['retest1SAMPLETEST210520_',num2str(xx),'_y_',num2str(yy),'_synthetic_350x350x5_den0p01_len20_cr0p7_NA_',num2str(NA),'_kerlt_',num2str(xyk),'_std_Gnoise_inc_',num2str(get_current),'.tif'])
close
end
toc
end
