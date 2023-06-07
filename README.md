# -----bdc-----
##Berkeley Distortion Calculator

This package calculates geometric optical distortion, by comparing observations to some reference frame, usually the PCU

This script is step 3 of 3.
Step 1) Take PCU observations with take_distortion_obs.py
Step 2) Extract starlist from these images with KAI
Step 3) Run bdc to generate a distortion solution from the starlists.


##------Running------
You can run the code in a terminal with:

python bdc_main.py config1.yml

Where config1.yml is the name of your configuration file.


##-------Required Files-----------

###.yml configuration file
This file contains the configuration settings for your run

###.lis Starlist files
Starlists generated by KAI from the .fits files

###.fits Images
The image files. Used for plotting

###.log file
The log file generated by take_distortion_obs.py. This is only used to get the Upper Z and Rotation positions. Once these are in the FITS header, it can be ignored.

### Flystar python package



##------YAML Configuration file:------

####nights : ['n1','n2']
The script can handle data from several different observing runs, with potentially different targets and parameters. These runs are called 'nights' in the script. The names in this list point to the sections later in the config file, and are used in naming outputs.

####instrument : 'OSIRSIS'
The instrument we are calculating a distortion model for. Currenly only handles 'OSIRIS'.

####output_dir      : '/u/mfreeman/work/PCU/script_testing/distortion_testing/outputs5/'
This directory will be generated, and outputs from bdc will be saved there.

####ref_instrument : 'PCU'	
The reference frame that we are aligning the observations to. Currently only handles 'PCU' and 'Hubble'.

####legendre_order : 5	
The order of the Legendre polynomial that forms the final distortion model.

####generate_new_fp : True 
If False, will attempt to load a previously generated four-parameter transformation. If True, will always generate a new one. Set to False to save time when re-running the script.		

####generate_reference_frame : False	
If False, will attempt to load a previously generated combined reference frame. If True, will always generate a new one. Set to False to save time when re-running the script.

####fp_iterations : 5 
Sets the number of times that the four-parameter transformation is repeated in a loop, to improve accuracy.

####ref_iterations : 5	
Number of times to repeat the combined reference frame generation.

####previous_distortion_dir : '/u/mfreeman/work/PCU/script_testing/distortion_testing/outputs3/'
The location of the previous distortion model for the instrument.

####previous_distortion_model : 'distortion_coefficients_r0_f3.txt'  #filename for the distortion coefficients. Must also have a file appended with 'inverse_' for the inverse transformation.
The file containing the previous distortion model for the instrument, as a list of legendre polynomial coefficients. There must also be an inverse model with 'inverse\_' prepended to the same filename.

####centred : True 		
If True, will add a translation to the final result to set the distortion at the centre of the frame to be zero. (This is an extra term, the Legendre polynomial is not changed) 

####use_flystar_velocity  : False 
If True, Flystar will be used to propogate the positions of the stars in the reference frame forward to the observation epoch using the velocities in the reference list. If False, the propogation is performed in this script.

####manually_average_star_positions : True  
When generating a new reference frame, the positions of stars that appear in multiple frames are averaged. If this parameter if False, the average calculated by Flystar is used. If True, the average is calculated in this script (testing has shown no difference)

####n1:
The name of the first 'night' of data. Must match an element from the parameter 'nights'.
####  target          : 'PCU'
There are some hard coded options for the target, with RA, Dec, and epoch. This should be deleted
####  ref_instrument  : 'PCU'  
The instrument which provided the reference frame. Currently only works with 'PCU'. Could eventually add 'Hubble' 
####  reference_file  : '/u/mfreeman/code/python/bdc/pinhole_grid_positions.txt' 
The reference frame starlist. Contains a starlist with X and Y positions
####  log_file        : '/u/mfreeman/work/PCU/script_testing/day3/PCU_20230413_060326.log'
The log file generated by take_distortion_obs.py. Only needed for the Upper Z and Rotation stage positions, once these are in the FITS header keywords the log file can be ignored.
####  starlist_dir    : '/u/mfreeman/work/PCU/script_testing/day3/raw/'      
The location of the starlist files generated by KAI
####  fits_dir        : '/u/mfreeman/work/PCU/script_testing/day3/raw/'      
The location of the .fits files (or \_flip.fits for OSIRIS)
####  bad_files       : []
List any bad starlist files here to skip them.
####  tform_dir       : '/u/mfreeman/work/PCU/script_testing/day3/tforms5/'  
This directory will be created, and the 4-parameter transformation guesses will be saved to here, and loaded from here on subsequent runs to save time.
####  mag_limits      : [-5,16]      
Stars outside this magnitude range are discarded. Also passed as limits to mosaic_to_ref()
####  mag_tolerance   : [3, 3, 3]
The matching magnitude tolerance passed to Flystar
####  target_RA        : 0         
When using an on-sky reference, the RA is calculated relative to this position.
####  target_DEC       : 0
####  single_fit      : True  
If True, will run Mosaic2Ref for each image individually. Only works for True
####  rad_tolerance   : [0.4, 0.3, 0.2]   
Matching radius tolerance passed to Flystar, in arcseconds.
####  slice           : [0,null]
If you want to run the distortion on a subset of the starlists, you can select them with a python slice object here. Use [0,null] to select all files.


##------Outputs------
Outputs will be saved in the directory in the configuration file. Outputs include:

#### Plots

#### dist_measures.txt
A list of all matched stars used to calculate the distortion model, for each iteration

#### distortion_coefficients.txt
The coefficients of the legendre polynomial transformation that forms the distortion model.

#### inverse_distortion_coefficients.txt
The inverse of the distortion model.

#### t_params.txt
The scale and rotation of the four-parameter transformation. Used to correct the plate scale and position angle when matching to on-sky data (not the PCU, do not have exact position angle)

#### iteration_residuals.txt
Statistics for the residual between the distortion-corrected stars, and the reference stars.



