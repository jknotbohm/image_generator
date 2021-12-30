# READ ME

*Repository for synthetic image generator mimicking confocal microscopy.*

This document is intended to guide the user in using the synthetic image generating algorithm.

Written by Mainak Sarkar under guidance of Jacob Notbohm, University of Wisconsin-Madison, 2020-2021. https://notbohm.ep.wisc.edu/

## Included files

### Main Scripts 

For image generation (general): To generate a sequence of synthetic images representing reference and current configurations of a 3D model of specimen (at focal plane) under any given
boundary condition, use `synth_confocal_imgen.m`

For image generation in artificial translation: To generate a sequence of images of a 3D specimen (at focal plane) under a prescribed rigid body translation along x and y, please 
use `synth_confocal_imgen_ATA.m`

### Scripts, data folders, function files and necessary MAT files

Hard-coded data of 3D PSF for any confocal microscope (to extract dimensionless xy intensity kernel at a given z): `intensity_kernel.mat` and `I_data_ND_Finest.m `

Mathematica script to generate the 3D PSF: `I_cal_ND_infinity_sum.nb` (This script uses  equation and non-dimensional parameters as suggested in the following reference: Born M, Wolf E. Principles of optics: electromagnetic theory of propagation, interference 
and diffraction of light. Elsevier; 2013 Jun 1.; Boas, David A., Constantinos Pitris, and Nimmi Ramanujam, eds. Handbook of biomedical optics. CRC press, 2016.)

`XY_inten_kernel.m` (this extracts 2D xy intensity kernel in dimensionless form from 3D PSF data)

`tiff_image_gen_rev_bground.m`

`drawLine_3D.m`

`conv_with_intensity.m`

`crop_network.m` (for rectangular images in general) 

`crop_network_modified.m` (a simplified version for square images): We recommend using this version as long as the synthetic image is square. 

The subfolder "3D_d0p025f20_Network_350x350x5_Seed_1003" contains a 3D numerical model of reference network (with and without inclusion) generated using a simulated annealing based algorithm 
(reference of algorithm: Grimmer, P., Notbohm, J.: Displacement propagation in fibrous networks due to local contraction. JBiomech Eng-T ASME 140(4) (2018))

Load-step-wise nodal displacement data (for network) are stored by ABAQUS in an ODB file during FE simulation. Those data were extracted to a MAT file; a sample is provided here 
to aid the user in running the image generator script: `sample_odb_file_name_network_DATA.mat`

There are 2 function files for step-wise reading of the nodal displacement data from the MAT file obtained post-extraction of an ODB file: `def_20data.m `
and `inc_20def_coords_rev.m`

`ata_def_nodes.m` (this function is instrumental in executing the prescribed rigid body translation of the 3D network)

### Sample confocal image of collagen

`sample_confocal_image.tif` This is a sample image obtained from real confocal microscopy, which was used to get an idea of the noise generated during the addition of artifical noise on synthetic images.

## Description 

### Using image generator algorithm with default image generating parameters

Please note the following default image generating parameters in our algorithm: NA = 0.7 
(chosen to be less than the objective's specified NA = 0.75 to account for  reduced resolution caused by larger than optimal pinhole), 20x air objective with n = 1, imaging in light of wavelength = 525 nm, camera sensor's 
pixel size = 6.5 µm. You can choose other values of these parameters by updating their input values in designated lines in `synth_imgen_confocal.m` and `synth_imgen_confocal_ATA.m`. 

Download the folder image_generator folder and set it as your current path. Also set this folder in the "set path" under home in task bar of MATLAB with sub-folders.

Following are step-by-step instructions to use this code.

#### For generating reference and a sequence of deformed images under a boundary condition
Step 1: Save the mat file `intensity_kernel.mat` in your current folder. 

OR,

Run the script `I_data_ND_Finest.m` once. This step is required in the very first run of the simulation only.

Step 2: Locate the MAT file that contains step-wise deformation data and load the same in `synth_imgen_confocal.m` script.

*Note:* To run the example considered in this document, please use the MAT file containing step-wise deformation data located here: 
`image_generator/imgen_function_files/sample_odb_file_name_network_DATA.mat` and save it in your current folder. By default this MAT file has been invoked in the main script 
`synth_imgen_confocal.m`.

Step 3: At the end of the script `synth_imgen_confocal.m`, write the name of tiff image files and the locations where you want them to be saved. A multi-tif image file is also generated, so make sure you name the ‘save path’ properly.

Step 4: Run `synth_imgen_confocal.m`. You can monitor the progress of image generation on the command prompt.

#### For reference and translated images in rigid body translation analysis
Step 1: Save the mat file `intensity_kernel.mat` in your current folder.

OR,

Run the script `I_data_ND_Finest.m` once. 

Step 2: In the script `synth_imgen_confocal_ATA.m`, specify x and y displacements.

Step 3: At the end of the script `synth_imgen_confocal.m` script, write the name of tiff image files and the locations where you want them to be saved. 

Step 4: Run `synth_imgen_confocal_ATA.m`.

### Example
To get an idea of how this works, please run any of the two main scripts (either general or the rigid body translation analysis) with all the default parameters. For the general one (`synth_confocal_imgen.m`), if you 
want the reference and final deformed state of the model, go to line number 30 and set get_current = [0 20] and in line number 194 set k = [0 get_current]. The given MAT file contains simulation results 
in 20 load steps: get_current = 0 corresponds to the reference and 20 corresponds to the final deformed state. You can customize get_current variable to get image of any intermediate configuration of the network.

### Sample synthetic images of network A, B, C, D
 
Sample synthetic images can be found here: https://uwmadison.box.com/s/i3tjcdkbe4kjg091vwhvi0faw4asi16a

