#this is the main script that calculates the distortion.

"""

Inputs:

Image starlists: We need to run starfinder on the images to get starlists. That could be included as a step here, or it could be a separate piece of code.
Either way, this script should take in the starlists.

Reference starlist: This will probably be the PCU grid, but could also be Hubble/Gaia catalogues. Make that an option.
Other target info? M15, M92, PCU grid

Initial transformation guesses: Very important for using the PCU. Maybe also important for Hubble/Gaia. I did end up using initial guesses. Might be hard to automate.

Create stacks: Combine images with the same PA into stacks with error measurements. Might need some intelligence here to pick the correct stacks.

Maybe include functions for f_tests etc to compare.

command looks like: 

python bdc_main.py -c config1.yml


To Do list:

Record residuals and other info. Save to text file? Or just print out at end.
Maybe make a dictionary of saved parameters?
I want the residual from each match I think.
Each combo of reference iteration and fourP iteration. I think I alreay saved these in a file.
Results dictionary with lists: four_p iteration, reference iteration, mean residual, weighted mean residual. Maybe I only need the residuals?
Maybe include number of matched stars too. Average distortion size. and standard deviation. Standard deviation of residuals.
Append all of these to lists.


Loading previous transformaiton guesses.
Need a proper system to do this. Should probably be a separate function?
For the PCU, just use the header params to calculate a guess each time. Should be good enough.
For a reference like Hubble, it is more complicated. Need to create a folder where the guesses can be stored.
And a clear naming convention so that files can be renamed. 
Also must work with no guess. Because the first pass through will not have any guesses
This happens again for the matched referenc frame. Although, I could possibly just use all the individual transformations?
This may not work, so save a list of transformations too.


Plots
Can really trim down. Only need ... some? Check the ones I used in my paper.
Difference vectors
Distortion model vectors
residuals.
Plate scale and rotations.


prepare_hubble_for_flystar needs to be flexible. Maybe rename everything?

I pass the globular cluster name in. This should be changed. I use it for the RA and Dec, then use 
The Hubble reference list has its RA and Dec transformed to be relative to the COO star. So I was getting the ra and Dec of the coo star.
I have the Ra and Dec of the coo star from GAIA, and the centre of the whole cluster from another paper.
I convert them relative to the centre of the cluster, using pixels relative to Coo star, and Coo star RA and Dec relative to cluster centre.
Plus a manual correction.
Solution ... 
This is a function in flystar because we want everything in coordinates relative to the target. That is not really necessary here. 
But we probably should still follow the convention.
I really need to look at that function, is it actually correct? I guess the reference is arbitrary, so it works either way.
The Hubble list x_0 is in the rectified Cartesian system wrt to the adopted center [arcsec]. So I think it is correct.
So that's why I don't need to change anything. I could probably just leave hubble as it is.
That will be inconsistent with my previous runs.
New plan:
Don't need to pass in Ra, Dec, or Target name. prepare_ref_for_flystar() doesn't need to translate hubble.
need prepare_PCU_for_flystar, which does something similar.
or just make the PCU list in the correct format anyway? It won't be in RA and Dec, and trying to convert it may be incorrect.
Do I need to make up the columns too? Magnitudes are a bit meaningless.
Maybe I could point to prepare_gaia for future changes.

"""
import argparse
import yaml
from astropy.table import Table
from astropy.io import ascii, fits
from flystar import analysis, align, transforms, starlists, plots
from kai.reduce import dar
from kai import instruments
import pickle
import os
import fnmatch
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib
# matplotlib.use('agg') #stops memory leak?
from matplotlib.colors import LogNorm, Normalize
import numpy as np
import math
import datetime
import copy
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", required = True, help = "Configuration file")
# parser.add_argument("-o", "--output", required = True, help = "Output folder")
args = parser.parse_args()
print("Config file is: {}".format(args.config))
with open(args.config, "r") as ymlfile:
	config = yaml.safe_load(ymlfile)

# print("Output folder is : {}".format(args.output))
print("Output folder is : {}".format(config['output_dir']))
plot_dir = config['output_dir'] + 'plots/'
os.makedirs(config['output_dir'],exist_ok=True)
os.makedirs(plot_dir,exist_ok=True)
refTable_current_filename = '{}refTable_current.pickle'.format(config['output_dir'])  #filename to save the refTable pickle
if config['instrument'] == 'OSIRIS':
	kai_instrument = instruments.OSIRIS()
	frame_size = [2048,2048] #could read out hdr['NAXIS1'] and hdr['NAXIS2']
else:
	print('Need to add a frame_size for {}'.format(config['instrument']))
frame_centre = [frame_size[0]/2,frame_size[1]/2]  
trim_not_used_in_trans=False
pcu_x_keyword = 'PCSFX'
pcu_y_keyword = 'PCSFY'
pcu_z_keyword = 'PCUZ'
pcu_r_keyword = 'PCUR'

def main():
	refTable_current = {}	#holds most recent refTable for each night #this is the reference table that keeps getting updated by the reference section.
	results = {}
	# old_distortion_model = None 
	distortion_model = distortion_from_coefficients(config['previous_distortion_dir']+config['previous_distortion_model']) #load old distortion model
	inverse_distortion_model = distortion_from_coefficients(config['previous_distortion_dir']+'inverse_'+config['previous_distortion_model'])
	if config['ref_iterations'] > 0:
		for ref_iteration in range(config['ref_iterations']):
			print('\n \n Ref Iteration {} \n'.format(ref_iteration))
			distortion_model, inverse_distortion_model, refTable_current = distortion_section(refTable_current, ref_iteration,results,distortion_model,inverse_distortion_model)
			new_refTable = reference_section(refTable_current,distortion_model,ref_iteration,results)
			refTable_current = new_refTable
	else:
		distortion_model,inverse_distortion_model, refTable_current = distortion_section(refTable_current,0,results,distortion_model,inverse_distortion_model)
		#outputs are printed to files.
	return distortion_model

def distortion_section(refTable_current,ref_iteration,results,current_distortion_correction=None,inverse_distortion_correction=None):
	#this section takes a reference table and some OSIRIS observations, matches them, and calculates a distortion model using Legendre polynomials.

	# current_distortion_correction = initial_distortion_correction
	# current_distortion_correction = None
	tab1_initial = {} #we only want to use the match table from the first four-parameter iteration
	for night in config['nights']:
		tab1_initial[night]=[]
		os.makedirs(config[night]['tform_dir'],exist_ok=True)


	# do I need the empty lists here in order to append to them?
	#if I am loading dist files then tab1_intial can break if a set of fp_iterations was incomplete - it won't run the first iteration to generate tab1_initial. Need to save as a file?
	for fp_iteration in range(config['fp_iterations']):

		matched_star_table_name = '{}dist_measures_r{}_f{}.txt'.format(config['output_dir'],ref_iteration,fp_iteration)
		matched_star_table = {}   #matched star lists from each image are all appended to this dictionary. Should be an astropy table? Lots of appending lists. Convert later?
		plate_scale_and_rotations = {}


		if config['generate_new_fp']:  		#if config['generate_new_fp']==True, we generate a new table of matched stars. If False, we load the previouly saved one.
			# obs_nights = [night for night in config]

			#-------------------------Night loop starts here-------------------------------
			for night in config['nights']:
				reference_data = fetch_reference(night)
				osiris_filenames = get_image_filenames(config[night]['starlist_dir'],config[night]['slice'])
				failed_images = []
				successful_images = []
				if night not in refTable_current.keys():  #if we have not generated a new reference frame, load the starting one.
					if config['use_flystar_velocity']:   #if True, use flystar to project positions. If False, manually project positions.
						refTable_H = prepare_ref_for_flystar(reference_data,config[night]['target_RA'],config[night]['target_DEC'],config[night]['target'],config[night]['ref_instrument'])
					else:
						starlist0 = load_osiris_file(config[night]['starlist_dir'] ,osiris_filenames[0]) #load 1 file to get the epoch.
						reference_data_p = project_pos(reference_data,starlist0,config[night]['ref_instrument'],config[night]['target'])
						refTable_H = prepare_ref_for_flystar(reference_data_p,config[night]['target_RA'],config[night]['target_DEC'],config[night]['target'],config[night]['ref_instrument'])
					
					refTable_current[night] = refTable_H.filled() #is this necessary?
					with open(refTable_current_filename, 'wb') as temp:
						pickle.dump(refTable_current, temp)                #sav

				for i, filename in enumerate(osiris_filenames):
					if filename in config[night]['bad_files']:
						print('{} {} flagged as bad, skipping'.format(i,filename))
						continue

					starlist = load_osiris_file(config[night]['starlist_dir'] ,filename)
					if config['instrument'] == 'OSIRIS':
						fitsfile = filename[:-10] + '.fits'     #for _stars.txt files
						#This fuction finds the fits file name by replacing _stars.txt with .fits
						#It assumes that the fits file names look like i230413_a003004_flip.fits and that the star_list file names look like i230413_a003004_flip_stars.txt
						#KAI output files may look like may look like ci211024_a005002_flip_0.8_stf.lis, use this:					
						#fitsfile = config[night]['fits_dir'] + filename[:-12] + '.fits'     #for _0.8_stf.lis files
					else:
						print('Need to add filename conventions for {}'.format(config['instrument']))

					PA = get_PA(config[night]['fits_dir']+fitsfile)
					# starlist = brightest_n(starlist,170)
					starlist = mag_cut(starlist,config[night]['mag_limits'][0],config[night]['mag_limits'][1])
					if not filename in config[night]['dont_trim']:
						starlist = edge_cut(starlist,5,frame_size)
					if len(starlist) == 0:					
						print(i,filename, '0 stars remaining after edge cut, skipping image')     #this needs to be outside the function.
						failed_images.append(filename[:-4])
						continue
					# refTable_t = trim_gaia(refTable_current[night],filename,PA)    
					#This function trims off stars that are outside the frame. Broken at the moment. Could still be useful? Will need different modes for PCU and on sky.
					refTable_t = refTable_current[night]
					refTable_tf = refTable_t.filled()		#unmask, required for quiver plots. Can this be removed?
					if config[night]['ref_instrument'] == 'PCU':
						refTable_d = refTable_tf[:]   #The PCU does not need to have DAR applied
					else:
						refTable_d = dar.applyDAR(config[night]['fits_dir']+fitsfile, refTable_tf, plot=False, instrument=osiris, plotdir=plot_dir + 'dar_d/'+ str(ref_iteration) + '/')

					try:

						transformation_guess_file = '{}tform_r{}_f{}_{}.pickle'.format(config[night]['tform_dir'],ref_iteration,fp_iteration,filename)

						if os.path.exists(transformation_guess_file):
							print('Loading transform guess')
							with open(transformation_guess_file, 'rb') as trans_file:
								transformation_guess = pickle.load(trans_file) 	
						else:
							if config[night]['ref_instrument'] == 'PCU':
								transformation_guess = calculate_PCU_guess(filename,night)							
							else:
								# transformation_guess = last_good_transform
								print('Not loading transform guess')
								transformation_guess = None

						# transformation_guess = calculate_PCU_guess(filename,night)							

						starlist_corrected = starlist[:]
						# if fp_iteration > 0:
						if (current_distortion_correction is not None) and (inverse_distortion_correction is not None):    #either this check, or the fp_iteration check.
							xt, yt = current_distortion_correction.evaluate(starlist_corrected['x'],starlist_corrected['y'])
							starlist_corrected['x'] = xt 
							starlist_corrected['y'] = yt 
						if config[night]['single_fit']:
							refTable_d_f = refTable_d #.filled(). not filling.
							msc = align.MosaicToRef(
								refTable_d_f, [starlist_corrected], iters=3,
							    dr_tol=config[night]['rad_tolerance'], dm_tol=config[night]['mag_tolerance'],
							    outlier_tol=[None, None, None],
								# trans_class=transforms.PolyTransform,
								trans_input=transformation_guess,
								trans_class=transforms.four_paramNW,
								trans_args=[{'order': 1}, {'order': 1}, {'order': 1}],
							    use_vel=False,
							    use_ref_new=False,
							    update_ref_orig=False, 
							    mag_trans=True,
							    mag_lim=config[night]['mag_limits'], #[6,16],
							    weights='both,std',
							    calc_trans_inverse=True,    
							    init_guess_mode='miracle', verbose=0) #6
							msc.fit()
							# tab1 = msc.ref_table
							tform = msc.trans_list
							tform_inv = msc.trans_list_inverse
							#I only want to save the transformation from the subsequent iterations, not tab1.
							# if fp_iteration == 0:
							# 	tab1_initial[night].append(msc.ref_table)
							# 	#save tab1_initial to a file. Once it has all nights in it?
							# tab1 = tab1_initial[night][i]
							
							j = 0
							tab1 = msc.ref_table
							if (current_distortion_correction is not None) and (inverse_distortion_correction is not None):
								tab1['x_orig'][:,j], tab1['y_orig'][:,j] = inverse_distortion_correction.evaluate(tab1['x_orig'][:,j], tab1['y_orig'][:,j])
	
							#
							
						else:
							j=i  # i=osiris index, j=gaia index
							print('Only set up for single_fit mode')
						ref_idx = np.where(tab1['ref_orig'] == True)[j]
						if len(ref_idx) <= 10:
							print(i, filename, 'Only', len(ref_idx),'matched, skipping')
							failed_images.append(filename[:-4])
							break  #return works too?
						with open(transformation_guess_file, 'wb') as temp:
							pickle.dump(tform, temp)

						print(i, filename, len(refTable_d), 'Reference stars,', len(starlist_corrected), 'OSIRIS stars,', len(ref_idx), 'matches')
						
						# px = tform_inv[j].px
						# py = tform_inv[j].py
						# theta = math.atan2(px[2],px[1])
						# scale = math.cos(theta) / px[1] #this may be a mistake, google says it's px[1]/cos(theta). Maybe that's why the inverse transformation works?
						
						px = tform[j].px
						py = tform[j].py
						theta = math.atan2(px[2],px[1])
						scale = px[1] / math.cos(theta)	

						plate_scale_and_rotations.setdefault('Filename',[]).append(filename)        #the setdefault() method adds to dictionary element, and creates it if it does not exits.
						plate_scale_and_rotations.setdefault('Scale',[]).append(scale)
						plate_scale_and_rotations.setdefault('Rotation',[]).append(math.degrees(theta))
						plate_scale_and_rotations.setdefault('PA',[]).append(PA)

						append_to_matched_star_table(matched_star_table,tab1,ref_idx,tform,tform_inv,night,i)

						x = tab1['x_orig'][ref_idx,j]
						y = tab1['y_orig'][ref_idx,j]
						# ids = tab1['name'][ref_idx]   #ID's of stars?
						Raref = tab1['x0'][ref_idx]
						Decref = tab1['y0'][ref_idx]						
						xref, yref = tform_inv[j].evaluate(Raref,Decref)						
						plot_scatter(starlist_corrected['x'],starlist_corrected['y'], fitsfile,'{}{}_scatter/'.format(plot_dir,night),night,fp_iteration,ref_iteration)
						plot_matched(x,y,xref,yref, fitsfile,'{}{}_matched/'.format(plot_dir,night),tab1['use_in_trans'][ref_idx],night,fp_iteration,ref_iteration)
						plot_quiver(x,y,xref,yref,fitsfile,'{}{}_quiver/'.format(plot_dir,night),tab1['use_in_trans'][ref_idx],night,fp_iteration,ref_iteration)
					
					except AssertionError as err:
						print(filename[:-4], 'Assertion error:')
						print(err)
						failed_images.append(filename[:-4])
						continue

					except ValueError as err:
						print(filename[:-4], 'Value error:')
						print(err)
						print('Usually mean no stars matched. Check dm_tol dr_tol')
						# print(err.__traceback__)
						failed_images.append(filename[:-4])
						raise
						continue

					successful_images.append(filename[:-4])

				#--------------------------------------------Night loop finishes here----------------------------------------------

				print('--------------------')
				print('Succeeded for', len(successful_images), 'files.')
				print('Failed for', len(failed_images), 'files:')
				print(failed_images)

				for key, value in plate_scale_and_rotations.items():
					plate_scale_and_rotations[key] = np.array(value)

				print(plate_scale_and_rotations)

				plate_scale_and_rotations['Difference'] = np.array(plate_scale_and_rotations['Rotation']) - np.array(plate_scale_and_rotations['PA'])
				transform_table = Table(plate_scale_and_rotations)
				ascii.write(transform_table,'{}t_params_{}_r{}.txt'.format(config['output_dir'],night,ref_iteration),format='fixed_width', overwrite=True) #may be inaccurate when using PCU
				# print('Mean scale =', np.mean(scales))
				# print('Mean rotation offset =', np.mean(plate_scale_and_rotations['Difference']))
				
				output_table = Table(matched_star_table) #elements are lists, not numpy arrays. Should be fine?
				ascii.write(output_table,matched_star_table_name,format='fixed_width', overwrite=True)
	

			#save tab1_initial to a file here? Now it has all nights in it. Can then always read in if trying to restart.
		else:
			print('b: Loading fit {}'.format(matched_star_table_name))
		
		distortion_data = ascii.read(matched_star_table_name,format='fixed_width')

		# use the parameters in output table to calculate a new distortion solution (check fit_legendre script)0
		# use that distortion solution to update correction_list
		if trim_not_used_in_trans:
			print('Trimming distortion_list to only include stars used in transformation')
			used_in_trans_B = np.array(distortion_data['UiT'].data) #copy of data, not pointer to
			used_in_trans_B = (used_in_trans_B == "True")  #Converts "True" strings to True booleans
			used_in_trans_B.tolist()
			distortion_data = distortion_data[used_in_trans_B]
			#I can probably compress this. Make the second line a list not a tuple?

		# flystar rejects outliers and sets their weights to zero, which sets their errors to inf. So we put them back here.
		ind_r = np.where(distortion_data['Weight'] == 0.0)[0] #index of stars rejected with outlier_rejection, have their weights set to 0.0
		distortion_data['Weight'][ind_r] = 1/np.sqrt((distortion_data['xe_IM'][ind_r]*0.01)**2 + (distortion_data['ye_IM'][ind_r]*0.01)**2 + distortion_data['xe_REF'][ind_r]**2 + distortion_data['ye_REF'][ind_r]**2)
		print('{} stars rejected due to Flystar outlier_rejection'.format(len(ind_r)))

		#maybe don't declare all these variables, just use the table, and update the weights in place.

		outliers, bad_names, m_distances = find_outliers(distortion_data['x_IM'],distortion_data['y_IM'],distortion_data['x_REF'],distortion_data['y_REF'],distortion_data['Ref_id'])
		include = ~outliers 
		print('Mahalanobis distance. Mean: {}, STD: {}'.format(np.mean(m_distances),np.std(m_distances)))

		print('Fitting Legendre Polynomial')
		legendre_transformation = transforms.LegTransform.derive_transform(distortion_data['x_IM'][include], 
																			distortion_data['y_IM'][include], 
																			distortion_data['x_REF'][include], 
																			distortion_data['y_REF'][include], 
																			config['legendre_order'], m=None, mref=None,init_gx=None, init_gy=None, weights=None, mag_trans=True
																			) # Defines a bivariate legendre tranformation from x,y -> xref,yref using Legnedre polynomials as the basis.
		x_coefficient_names = []
		x_coefficient_values = []
		y_coefficient_names = []
		y_coefficient_values = []
		for param in legendre_transformation.px.param_names:
			# print(getattr(legendre_transformation.px, param))
			a = getattr(legendre_transformation.px, param)
			x_coefficient_names.append(a.name)
			x_coefficient_values.append(a.value)
		for param in legendre_transformation.py.param_names:
			a = getattr(legendre_transformation.py, param)
			y_coefficient_names.append(a.name)
			y_coefficient_values.append(a.value)
		output_table = Table([x_coefficient_names,x_coefficient_values, y_coefficient_names,y_coefficient_values], names = ('px_name','px_val','py_name','py_val'),)
		ascii.write(output_table,'{}distortion_coefficients_r{}_f{}.txt'.format(config['output_dir'],ref_iteration,fp_iteration),format='fixed_width', overwrite=True)


		inverse_legendre_transformation = transforms.LegTransform.derive_transform(distortion_data['x_REF'][include], distortion_data['y_REF'][include], distortion_data['x_IM'][include], distortion_data['y_IM'][include], 
																			config['legendre_order'], m=None, mref=None,init_gx=None, init_gy=None, weights=None, mag_trans=True
																			) #The inverse of the distortion correction: adds distortion to an undistorted image.
		x_coefficient_names = []
		x_coefficient_values = []
		y_coefficient_names = []
		y_coefficient_values = []
		for param in inverse_legendre_transformation.px.param_names:
			a = getattr(inverse_legendre_transformation.px, param)
			x_coefficient_names.append(a.name)
			x_coefficient_values.append(a.value)
		for param in inverse_legendre_transformation.py.param_names:
			a = getattr(inverse_legendre_transformation.py, param)
			y_coefficient_names.append(a.name)
			y_coefficient_values.append(a.value)
		output_table = Table([x_coefficient_names,x_coefficient_values, y_coefficient_names,y_coefficient_values], names = ('px_name','px_val','py_name','py_val'),)
		ascii.write(output_table,'{}inverse_distortion_coefficients_r{}_f{}.txt'.format(config['output_dir'],ref_iteration,fp_iteration),format='fixed_width', overwrite=True)

		distortion_plot_print_save(distortion_data,include,legendre_transformation,inverse_legendre_transformation,fp_iteration,ref_iteration,results)


	return legendre_transformation, inverse_legendre_transformation, refTable_current #returns the final distortion correction




def reference_section(refTable_current,distortion_model,ref_iteration,results):
	#this function applies a distortion solution to the OSIRIS lists, and combines them into a single ref list to return (per night)
	if len(refTable_current) == 0:
		with open(refTable_current_filename, 'rb') as temp:
			refTable_current = pickle.load(temp)   #refTable_current should be populated from the distoriton_section. Can also read it in as a file here, maybe for doing the final loop.

	#------------------------Night loop 1 starts here---------------------------

	# obs_nights = [night for night in config]
	# print(obs_nights)
	for night in config['nights']:

		osiris_filenames = get_image_filenames(config[night]['starlist_dir'],config[night]['slice'])
		print(len(osiris_filenames), 'OSIRIS images')

		#perhaps make this next section possibly generate a new reference frame, then always load from file. For consistency.
		#make a function to choose the ref file? It's just checking if there's a previous one, can probably do one step.

		# combined_ref_filename = config['output_dir'] + 'combined_ref_table_' + str(ref_iteration) + '.txt'
		combined_ref_filename = '{}combined_ref_table_{}_r{}.pickle'.format(config['output_dir'],night,ref_iteration)
		failed_images = []
		successful_images = []

		if config['generate_reference_frame'] or (not os.path.exists(combined_ref_filename)):
			list_of_starlists = []
			print('Generating new combined refTable')

			for i, filename in enumerate(osiris_filenames):
				if filename in config[night]['bad_files']:
					print('{} {} flagged as bad, skipping'.format(i,filename))
					continue
					# failed_images.append(filename[:-4])
				print('{} {} applying distortion correction'.format(i,filename))
				
				if config['instrument'] == 'OSIRIS':
					fitsfile = filename[:-10] + '.fits'     #for _stars.txt files
					#This fuction finds the fits file name by replacing _stars.txt with .fits
					#It assumes that the fits file names look like i230413_a003004_flip.fits and that the star_list file names look like i230413_a003004_flip_stars.txt
					#KAI output files may look like may look like ci211024_a005002_flip_0.8_stf.lis, use this:					
					#fitsfile = config[night]['fits_dir'] + filename[:-12] + '.fits'     #for _0.8_stf.lis files
				else:
					print('Need to add filename conventions for {}'.format(config['instrument']))

				# fitsfile = config[night]['fits_dir'] + filename[:-12] + '.fits'
				PA = get_PA(config[night]['fits_dir']+fitsfile)
				starlist = load_osiris_file(config[night]['starlist_dir'] ,filename)
				# starlist = brightest_n(starlist,170)
				starlist = mag_cut(starlist,config[night]['mag_limits'][0],config[night]['mag_limits'][1])
				if not filename in config[night]['dont_trim']:
					starlist = edge_cut(starlist,5,frame_size)
				if len(starlist) == 0:
					print(i,filename, '0 stars remaining after edge cut, skipping image')
					failed_images.append(filename[:-4])
					continue

				#load_and_prepare_data does something similar, but also trims the reference. Maybe split that up and use the function.
				#or just don't use the function, it's only a few lines.

				xt, yt = distortion_model.evaluate(starlist['x'],starlist['y'])    #apply distortion correction

				if config['centred']:
					xc, yc = distortion_model.evaluate(frame_centre[0],frame_centre[1])
					xc -= frame_centre[0]
					yc -= frame_centre[1]						
				else:
					xc = 0
					yc = 0
				starlist['x'] = xt - xc
				starlist['y'] = yt - yc

				# plt.figure(num=4,figsize=(6,6),clear=True)			
				if config[night]['ref_instrument'] != 'PCU':
					starlist = dar.removeDAR(config[night]['fits_dir']+fitsfile,starlist, plot=False, instrument=osiris, plotdir=plot_dir + 'dar_r/'+ str(ref_iteration) + '/')

				list_of_starlists.append(starlist)

				# successful_images.append(filename[:-4])

			print('--------------------')
			print('Completed applying distortion corrections')

			# The transformation for each starlist is saved to this file, and can be loaded as an initial guess.
			reference_transformation_guess_file = '{}reference_section/tform_{}_r{}.pickle'.format(config[night]['tform_dir'],night,ref_iteration)
			os.makedirs('{}reference_section/'.format(config[night]['tform_dir']),exist_ok=True)

			# If we don't have the above file, we can use all the individual transformations from the distortion section as a guess.
			transformation_guess_file_list = ['{}tform_r{}_f{}_{}.pickle'.format(config[night]['tform_dir'],ref_iteration,int(config['fp_iterations']-1),i) for i in osiris_filenames]
			# transformation_guess_file = '{}tform_{}_{}_{}.pickle'.format(config[night]['tform_dir'],ref_iteration,fp_iteration,filename)

			use_individual_trans_guesses = False
			if os.path.exists(reference_transformation_guess_file):
				print('Loading transform from {}'.format(reference_transformation_guess_file))
				with open(reference_transformation_guess_file, 'rb') as trans_file:
					trans_list_temp = pickle.load(trans_file)

				if len(trans_list_temp) == len(list_of_starlists):
					print('Transform guess has correct number of lists {}'.format(len(trans_list_temp)))
					reference_transformation_guess = trans_list_temp
				else:
					use_individual_trans_guesses = True
			else:
				use_individual_trans_guesses = True

			if use_individual_trans_guesses:
				print('Loading old transform')
				if  all([os.path.exists(f) for f in transformation_guess_file_list]):
					reference_transformation_guess = []
					for i, tform_filename in enumerate(transformation_guess_file_list):
						with open(tform_filename, 'rb') as trans_file:
							reference_transformation_guess.extend(pickle.load(trans_file))					#initial guess for the transformation
																						#each file is a list with one element, so using .extend
				else:
					# trans_list = last_good_transform
					print('No guesses found for reference transformation')
					reference_transformation_guess = None


			#Do I need to be able to pass these parameters in, or should they stay fixed? My testing found False/False to work best.
			set_use_ref_new = False
			set_update_ref_orig = False
			dr_tolerance = [0.5,0.5,0.5] #these should probably be config parameters
			dm_tolerance = [2,2,2]

			print('Creating combined reference frame with MosaicToRef')
			msc = align.MosaicToRef(refTable_current[night],
				list_of_starlists, iters=3,
				dr_tol=dr_tolerance, dm_tol=dm_tolerance,
				outlier_tol=[None,None,None],
				trans_input=reference_transformation_guess,
				trans_class=transforms.four_paramNW,
				trans_args=[{'order': 1}, {'order': 1}, {'order': 1}],
				use_vel=config['use_flystar_velocity'],
				mag_trans=True,
				mag_lim=config[night]['mag_limits'],
				weights=None,
				use_ref_new = set_use_ref_new,
				update_ref_orig = set_update_ref_orig,
				calc_trans_inverse=True,    
				init_guess_mode='miracle', verbose=0)


			msc.fit()
			# tab2 = msc.ref_table
			with open(reference_transformation_guess_file, 'wb') as temp:
				pickle.dump(msc.trans_list, temp)
			print('Completed reference table matching')
			# refTable_a = msc.ref_table
			with open(combined_ref_filename, 'wb') as temp:
				pickle.dump(msc.ref_table, temp)

		else:
			print('reference_section: Loading previous combined refTable')

		#----------- end creation of combined_ref_file -------------------

		with open(combined_ref_filename, 'rb') as combined_ref_file:
			refTable_a = pickle.load(combined_ref_file)
		

		#average the star points manually here? I think they are all saved in columns of refTable_a
		


		refTable_b = Table([refTable_a['name']], names=['name'], masked=False)
		if config['manually_average_star_positions']:
			refTable_b['x0'] = np.nanmean(refTable_a['x'],axis=1)
			refTable_b['y0'] = np.nanmean(refTable_a['y'],axis=1)
		else:
			refTable_b['x0'] = refTable_a['x0']
			refTable_b['y0'] = refTable_a['y0']
		refTable_b['x0e'] = refTable_a['x0e']
		refTable_b['y0e'] = refTable_a['y0e']
		refTable_b['vx'] = refTable_a['vx']
		refTable_b['vy'] = refTable_a['vy']
		refTable_b['vxe'] = refTable_a['vxe']
		refTable_b['vye'] = refTable_a['vye']
		refTable_b['t0'] = refTable_a['t0']
		refTable_b['m0'] = refTable_a['m0']
		refTable_b['m0e'] = refTable_a['m0e']
		# refTable_b['use_in_trans'] = refTable_a['use_in_trans']		#this causes a problem. If this column is present, align.py flags stars that are not used as not 'Keepers', and their weight is left as zero.
																		#Potential change: copy the whole table and delete the 'use_in_trans' column.

		if trim_not_used_in_trans:
			print('Trimming refTable_b to only include stars used in transformation')
			refTable_b = refTable_b[refTable_a['use_in_trans']]			#using refTable_a due to the problem above with refTable_b

		temp, idx1, idx2 = np.intersect1d(refTable_b['name'].astype('str'),refTable_current[night]['name'],return_indices=True)

		reference_section_plots(refTable_b,refTable_current)

		residuals_r = np.hypot(refTable_current[night]['x0'][idx2]-refTable_b['x0'][idx1],refTable_current[night]['y0'][idx2]-refTable_b['y0'][idx1])
		results.setdefault('mean_residuals_r',[]).append(np.mean(residuals_r))
		print('ref_iteration {} mean residuals_r = {}'.format(ref_iteration,results['mean_residuals_r']))

		refTable_current[night] = Table(refTable_b,copy=True)

	return refTable_current








def find_outliers(x,y,xref,yref,star_id,verbose=False):
	print('Finding outliers')
	#scan grid over field. Flag any points as outliers.  Service(2016) did 205x205 pixel bins. removing 3 sigma outliers.
	dx = x - xref
	dy = y - yref
	# bins = np.arange(0,2048+1,256)
	xbins = np.arange(0,frame_size[0]+1,int(frame_size[0]/8))
	ybins = np.arange(0,frame_size[1]+1,int(frame_size[1]/8))
	x_binned = np.digitize(x,xbins)
	y_binned = np.digitize(y,ybins)

	outflag = np.zeros(len(x),dtype=bool)
	bad_stars = []
	m_distances = []

	looping = True

	if looping:
		previous_sum = -1
		loopcount = -1
		while sum(outflag) > previous_sum:	
			previous_sum = sum(outflag)
			loopcount += 1
			if loopcount >= 20:
				if verbose:
					print("Aborting after 20 loops")
				break		
			for i, xb in enumerate(xbins):
				for j, yb in enumerate(ybins):
					
					if i == 0 or j == 0:
						continue

					include = ~outflag
					idx = np.nonzero((x_binned[include] == i) & (y_binned[include] == j))[0]
					
					if len(dx[include][idx]) < 3:
						if verbose:
							print('{} stars in box, skipping'.format(len(dx[include][idx])))
						continue

					distances = mahalanobis(dx[include][idx],dy[include][idx])
					m_distances.extend(distances)
					out = np.nonzero(distances >= 3)[0]
					# outflag[include][idx[out]] = True  #boolean   #this doesn't work becasue we get a copt of outflag, not a view.
					outflag[np.where(include)[0][idx[out]]] = True				
					bad_stars.extend(star_id[include][idx[out]])

			if verbose:
				print('Loop', loopcount, ', Outliers removed =', sum(outflag))

	else:
		for i, xb in enumerate(xbins):
			for j, yb in enumerate(ybins):
				
				if i == 0 or j == 0:
					continue
				
				include = ~outflag
					# print(idx)
				idx = np.nonzero((x_binned == i) & (y_binned == j))[0]


				distances = mahalanobis(dx[idx],dy[idx])
				m_distances.extend(distances)
				out = np.nonzero(distances >= 2)[0]
				outflag[idx[out]] = 1  #boolean
				print(outflag[idx[out]])			
				bad_stars.extend(star_id[idx[out]])

				# else:
				# 	stdx = np.std(dx[idx])
				# 	stdy = np.std(dy[idx])
				# 	stdr = np.sqrt(stdx**2 + stdy**2)

				# 	meanx = np.mean(dx[idx])
				# 	meany = np.mean(dy[idx])

				# 	r = np.sqrt((dx[idx]-meanx)**2 + (dy[idx]-meany)**2)

				# 	rx = abs(dx[idx]-meanx)
				# 	ry = abs(dy[idx]-meany)
				# 	out = np.nonzero((rx > 3*stdx) | (ry > 3*stdy))[0]
	print('Outliers found')
	# bad_names.sort()    #probably don't need this section on bad reference stars.
	# bad_count = Counter(bad_names)
	# all_count = Counter(star_id)
	# bad_stars = [k for (k,v) in bad_count.items() if v > 11]
	return outflag, bad_stars, m_distances	

def mahalanobis(dx,dy):
	mu = np.array([np.mean(dx),np.mean(dy)])
	x = np.column_stack((dx,dy))
	S = np.cov(x.T)
	SI = np.linalg.inv(S)
	D = np.diag(np.sqrt(np.dot(np.dot((x-mu), SI), (x-mu).T)))
	return D





def fetch_reference(night):				#loads hubble data from hst1pass
	# filelist = os.listdir(nightDir)
	# hst1pass_file = fnmatch.filter(filelist,'i*.xymrduvqpk')[0]
	print('Loading', config[night]['reference_file'])
	# hubble_table = ascii.read(nightDir + hst1pass_file,format='commented_header',guess=False,delimiter='\s',comment='#',header_start=166,data_start=168,) #lines start at 0
	# hubble_table = ascii.read(nightDir + hst1pass_file,format='basic',names=['x','y','m','r','d','u','v','q','p','k'])
	if config[night]['ref_instrument'] == 'Hubble':
		column_names = ['r', 'x_0', 'y_0', 'pmx', 'pmy', '1sig_pmx', '1sig_pmy', 'x_M', 'y_M', 'Delta_x', 'Delta_y', 
						'1sig_pmx_mod', '1sig_pmy_mod', 'm_F606W', 'm_F814W', 'rms_F606W', 'rms_F814W', 'qfit_F606W', 'qfit_F814W', 
						'red_chi2_pmx', 'red_chi2_pmy', '1sig_intrcp_pmx', '1sig_intrcp_pmy', 'time_bsln', '1sig_intrcp_pmx_mod', '1sig_intrcp_pmy_mod', 
						'pm_ref', 'N_found', 'N_used', 'ID', 'delta_pmx', 'delta_pmy']
		ref_table = ascii.read(config[night]['reference_file'],format='basic',names=column_names)
	
	elif config[night]['ref_instrument'] == 'PCU':
		ref_table = ascii.read(config[night]['reference_file'],format='fixed_width')
	else:
		print("ref_instrument {} not supported for fetch_reference. Must be 'PCU' or 'Hubble'".format(config[night]['ref_instrument']))
		sys.exit(0)

	return ref_table

def prepare_ref_for_flystar(ref, ra, dec, target, instrument, targets_dict=None, match_dr_max=0.2):   
	#copied from prepare_gaia_for_flystar(). #year=2006 by default?
	"""
	Take a Gaia table (from astroquery) and produce a new table with a tangential projection
	and shift such that the origin is centered on the target of interest. 
	Convert everything into arcseconds and name columns such that they are 
	ready for FlyStar input.
	"""


	# Master frame in UVIS pixel scale (40 mas/pix). 
	# Cluster's center (from Goldsbury et al. 2010, AJ, 140, 1830, as updated 
	#  in 2011) (RA, Dec)_J000.0 = (21:29:58.33, +12:10:01.20) 
	#  Cluster's center on the master frame: (4989.28, 5017.50) [pix]
	# ra and dec are in decimal degrees. Shift so that they are at zero 
	
	if instrument == 'PCU':

		arbitrary_pinhole_magnitude = 9.0
		front_focus_scale = 1.385 #milliarcseconds per micron, which equals arcseconds per mm.  
		#Measurements with pcu give 138.5 pixels per mm -> 1385 mas per mm -> 1.385 arcsec per mm.


		ref_new = Table([ref['id']], names=['name'], masked=False)

		ref_new['x0'] = (ref['x'] * front_focus_scale) *-1            #In arcseconds. *-1 to get East positive
		ref_new['y0'] = ref['y'] * front_focus_scale

		ref_new['x0e'] = 0.001 * front_focus_scale * np.ones(len(ref_new))   #Advance reproductions states pinhole tolerance is +- 1 micron.
		ref_new['y0e'] = 0.001 * front_focus_scale * np.ones(len(ref_new))

		ref_new['vx'] = np.zeros(len(ref_new))    #arcsec per year
		ref_new['vy'] = np.zeros(len(ref_new))
		ref_new['vxe'] = np.zeros(len(ref_new))     #arcsec per year
		ref_new['vye'] = np.zeros(len(ref_new))     

		ref_new['t0'] = 2023 * np.ones(len(ref_new))
		ref_new['m0'] = arbitrary_pinhole_magnitude * np.ones(len(ref_new))    #actual pinhole brightness should be equal after flatfielding. Set arbitrary value of 9 for calculations.  

		return ref_new

	elif instrument == 'Hubble':

		if target == 'm15':
			cluster_centre = SkyCoord('21:29:58.33', '12:10:01.20', unit=(u.hourangle, u.deg), frame='icrs')     #hard coded m15
			ra_centre = cluster_centre.ra.degree     # in decimal degrees
			dec_centre = cluster_centre.dec.degree   # in decimal degrees
			ra_diff = (ra_centre-ra)*3600     #degrees to arcseconds
			dec_diff = (dec_centre-dec)*3600   

			manual_ra =  -(13.19) * osiris_scale    #manual offset to align images.
			manual_dec = -(15.15) * osiris_scale
			ra_diff += manual_ra
			dec_diff += manual_dec
		elif target == 'm92':
			cluster_centre = SkyCoord('17:17:07.39', '43:08:09.40', unit=(u.hourangle, u.deg), frame='icrs')     #hard coded m92
			ra_centre = cluster_centre.ra.degree     # in decimal degrees
			dec_centre = cluster_centre.dec.degree   # in decimal degrees
			ra_diff = (ra-ra_centre)*3600     #degrees to arcseconds
			dec_diff = (dec-dec_centre)*3600   

			manual_ra =  -(0) * osiris_scale    #manual offset to align images.
			manual_dec = -(0) * osiris_scale
			ra_diff += manual_ra
			dec_diff += manual_dec
		else:
			ra_diff = ra*3600     #degrees to arcseconds
			dec_diff = dec*3600   

			# print('{} can only be be m15 or m92'.format(target))

		# cos_dec = np.cos(np.radians(dec))
		# x = (ref['r'] - ra) * cos_dec * 3600.0   # arcsec
		# y = (ref['d'] - dec) * 3600.0           # arcsec
		#xe = ref['ra_error'] * cos_dec / 1e3      # arcsec
		#ye = ref['dec_error'] / 1e3               # arcsec

		# need to make a column of names. Maybe source ID too? Or can I skip it.

		# ref['source_id'] = range(len(ref['x_0']))

		# ref_new = table.Table([ref['source_id'].data.astype('S19')], names=['name'], masked=False)
		ref_new = Table([ref['ID']], names=['name'], masked=False)

		# ref_new['x0'] = x * -1.0
		# ref_new['y0'] = y
		#ref_new['x0e'] = xe
		#ref_new['y0e'] = ye
		ref_new['x0'] = (ref['x_0']-ra_diff) *-1            #*-1 to get East positive
		ref_new['y0'] = ref['y_0']-dec_diff
		ref_new['x0e'] = ref['1sig_intrcp_pmx'] * 40/1000   #40mas/pix to arcseconds
		ref_new['y0e'] = ref['1sig_intrcp_pmy'] * 40/1000 

		# Also convert the velocities. Note that gaia PM are already * cos(dec)
		#ref_new['vx'] = ref['pmra'].data * -1.0 / 1e3 # asec/yr
		#ref_new['vy'] = ref['pmdec'].data / 1e3
		#ref_new['vxe'] = ref['pmra_error'].data / 1e3
		#ref_new['vye'] = ref['pmdec_error'].data / 1e3
		ref_new['vx'] = ref['pmx'] / 1000    #pmx is in mas per year. Converting to arcsec per year
		ref_new['vy'] = ref['pmy'] / 1000
		ref_new['vxe'] = ref['1sig_pmx'] / 1000     #arcsec per year
		ref_new['vye'] = ref['1sig_pmy'] / 1000   

		ref_new['t0'] = ref['time_bsln']  #is this correct? time_bsln seems to be 6 years. Set use_flystar_velocity = False for now.
		# ref_new['source_id'] = ref['ID']

		# Find sources without velocities and fix them up.
		# idx = np.where(ref['pmy'].mask == True)[0]
		# ref_new['vx'][idx] = 0.0
		# ref_new['vy'][idx] = 0.0
		# ref_new['vxe'][idx] = 0.0
		# ref_new['vye'][idx] = 0.0

		ref_new['m0'] = ref['m_F814W']
		# ref_new['me'] = 1.09/ref['phot_g_mean_flux_over_error']
		# ref_new['parallax'] = ref['parallax']
		# ref_new['parallax_error'] = ref['parallax_error']

		# Set the velocities (and uncertainties) to zero if they aren't measured.
		idx = np.where(np.isnan(ref_new['vx']) == True)[0]
		ref_new['vx'][idx] = 0.0
		ref_new['vxe'][idx] = 0.0
		ref_new['vy'][idx] = 0.0
		ref_new['vye'][idx] = 0.0

		ref_new = ref_new.filled()  #convert masked colunms to regular columns


		if targets_dict != None:    #This section renames stars in ref_new using names from targets_dict.
			for targ_name, targ_coo in targets_dict.items():
				dx = ref_new['x0'] - (targ_coo[0] * -1.0)
				dy = ref_new['y0'] - targ_coo[1]
				dr = np.hypot(dx, dy)

				idx = dr.argmin()

				if dr[idx] < match_dr_max:
					ref_new['name'][idx] = targ_name
					print('Found match for: ', targ_name)

		return ref_new


	else:
		print("Instrument {} not recognized, must be 'PCU' or 'Hubble'".format(instrument))
		sys.exit(0)


def get_image_filenames(directory,file_slice):
	filelist = os.listdir(directory)
	# lis_files = fnmatch.filter(filelist,'ci*.lis')
	lis_files = fnmatch.filter(filelist,'i*_stars.txt')  #these may need to be converted to .lis files?
	lis_files.sort()
	sl = slice(file_slice[0],file_slice[1])
	lis_files = lis_files[sl]	
	return lis_files

def project_pos(refData_o, starlist,ref_instrument,target):
	#Assume that all osiris images have the same epoch.
	#GAIA['ra'] and GAIA['dec'] are in degrees, GAIA['pmra'] and GAIA['pmdec'] are in mas/year
	#hubble['x_0'] is in arcseconds, hubble['pmx'] is in mas per year
	refData = refData_o[:]   #make a copy of the original array
	osirisT = starlist['t'][0]
	if ref_instrument == 'Gaia':
		dt = osirisT - refData['ref_epoch']
		print('Projecting Gaia {:.3f} years'.format(dt[0]))
		idra=~refData['pmra'].mask    #select not masked
		iddec=~refData['pmdec'].mask 	
		refData['ra'][idra] = refData['ra'][idra] + refData['pmra'][idra]/(3600*1000) * dt[idra]        #mas to degrees
		refData['dec'][iddec] = refData['dec'][iddec] + refData['pmdec'][iddec]/(3600*1000) * dt[iddec]
		#need to include error, like hubble below.
		refData['ref_epoch'] = osirisT  #??

	elif ref_instrument =='Hubble':
		if target == 'm15':
			#m15 reference time is J2006.334577
			hubbleT = 2006.334577 								#hard coded m15
		elif target == 'm92':
			#m92 reference time
			hubbleT = 2006.277936
		# print(osirisT)	
		dt = osirisT - hubbleT
		print('Projecting Hubble {:.3f} years, {} to {}'.format(dt,hubbleT,osirisT))
		refData['x_0'] = refData['x_0'] + refData['pmx']*dt/1000      #mas to arcsec
		refData['y_0'] = refData['y_0'] + refData['pmy']*dt/1000
		#proper motion error is in mas per year. convert to pixels by dividing by 40mas/pix
		refData['1sig_intrcp_pmx'] = np.hypot(refData['1sig_intrcp_pmx'],  refData['1sig_pmx']/40*dt)
		refData['1sig_intrcp_pmy'] = np.hypot(refData['1sig_intrcp_pmy'],  refData['1sig_pmy']/40*dt)
		refData['ref_epoch'] = osirisT  # update the epoch to match osiris
	elif ref_instrument == 'PCU':
		print("Don't need to project PCU reference epoch")
	else:	 
		print('Need projection calculation for {}'.format(ref_instrument))
	return refData	


# def load_transformation_guess():
# 	if reference_instrument == 'Hubble':
# 		tform_file_1 = '{}tform_{}_{}_{}.p'.format(config[night]['tform_dir'],ref_iteration,fp_iteration,filename)	#don't need to enter the night, because the filename includes the date and so is unique
# 		tform_file_2 = '{}tform_{}_{}.p'.format(transform_files_location,ref_iteration,filename)
# 		tform_file_3 = '/u/mfreeman/work/d/transform_files/hubble/tform_{}.p'.format(filename)

# 	elif reference_instrument == 'GAIA':
# 		#out of date
# 		tform_file_1 = '{}tform_{}_{}.p'.format(config[night]['tform_dir'],ref_iteration,filename)
# 		tform_file_2 = '{}tform_{}_{}.p'.format(transform_files_location,fitmode,ref_iteration,filename)
# 		tform_file_3 = '/u/mfreeman/work/d/transform_files/gaia/tform_{}.p'.format(filename)

# 	if os.path.exists(tform_file_1):
# 		print('Loading correct transform guess')
# 		with open(tform_file_1, 'rb') as trans_file:
# 			trans_list = pickle.load(trans_file) 					#initial guess for the transformation
# 	elif os.path.exists(tform_file_2):
# 		print('Loading default transform guess')
# 		with open(tform_file_2, 'rb') as trans_file:
# 			trans_list = pickle.load(trans_file)
# 	elif os.path.exists(tform_file_3):
# 		print('Loading old transform guess')
# 		with open(tform_file_3, 'rb') as trans_file:
# 			trans_list = pickle.load(trans_file)
# 	else:
# 		# trans_list = last_good_transform
# 		print('Not loading transform guess')
# 		trans_list = None	
# 	return trans_list


def append_to_matched_star_table(matched_star_table,tab1,ref_idx,tform,tform_inv,night,i):
	j=0 #?? does this ever vary?

	# ids1 = np.where(tab1['name_in_list'] == 's1')[0]
	Ox = tab1['x_orig'][ref_idx,j]
	Oy = tab1['y_orig'][ref_idx,j]
	GRa = tab1['x0'][ref_idx]
	GDec = tab1['y0'][ref_idx]
	ids = tab1['name'][ref_idx]   #ID's of stars?
	Gx, Gy = tform_inv[j].evaluate(GRa,GDec)
	ORa, ODec = tform[j].evaluate(Ox,Oy)
	w = tab1['w'][ref_idx]
	print('Number of times weight=0.0:', (w==0.0).sum())
	# print(np.median(w))
	# matched_star_table['x_O'].extend(Ox)
	# matched_star_table['y_O'].extend(Oy)
	# matched_star_table['x_G'].extend(Gx)
	# matched_star_table['y_G'].extend(Gy)
	# matched_star_table['weights'].extend(w)
	# matched_star_table['idnums'].extend(ids)
	# matched_star_table['starlist_num'].extend([i] * len(ids))
	# matched_star_table['night_col'].extend([night] * len(ids))
	# matched_star_table['xe_O'].extend(tab1['xe_orig'][ref_idx,j])
	# matched_star_table['ye_O'].extend(tab1['ye_orig'][ref_idx,j])
	# matched_star_table['xe_G'].extend(tab1['x0e'][ref_idx])
	# matched_star_table['ye_G'].extend(tab1['y0e'][ref_idx])
	# matched_star_table['used_in_trans_1'].extend(tab1['use_in_trans'][ref_idx])
	# matched_star_table['Ra_O'].extend(ORa)
	# matched_star_table['Dec_O'].extend(ODec)
	# matched_star_table['Ra_G'].extend(GRa)
	# matched_star_table['Dec_G'].extend(GDec)

	# setdefault: if the key exists, return the value. Otherwise, create the key, set value to [], then return it. Either case can then be extended.
	matched_star_table.setdefault('x_IM',[]).extend(Ox) 		#image star x position (pixels)
	matched_star_table.setdefault('y_IM',[]).extend(Oy)
	matched_star_table.setdefault('x_REF',[]).extend(Gx) 		#reference star x position (pixels) (calculated using inverse transformation)
	matched_star_table.setdefault('y_REF',[]).extend(Gy)
	matched_star_table.setdefault('Weight',[]).extend(w)
	matched_star_table.setdefault('Ref_id',[]).extend(ids)
	matched_star_table.setdefault('Frame_num',[]).extend([i] * len(ids))
	matched_star_table.setdefault('Night',[]).extend([night] * len(ids))
	matched_star_table.setdefault('xe_IM',[]).extend(tab1['xe_orig'][ref_idx,j]) 	#image star x error (pixels)
	matched_star_table.setdefault('ye_IM',[]).extend(tab1['ye_orig'][ref_idx,j]) 
	matched_star_table.setdefault('xe_REF',[]).extend(tab1['x0e'][ref_idx]) 		#reference star x error (pixels)
	matched_star_table.setdefault('ye_REF',[]).extend(tab1['y0e'][ref_idx])
	matched_star_table.setdefault('Used_in_trans',[]).extend(tab1['use_in_trans'][ref_idx])
	matched_star_table.setdefault('RA_IM',[]).extend(ORa) 			#image star RA (calculated using transformation)
	matched_star_table.setdefault('Dec_IM',[]).extend(ODec)
	matched_star_table.setdefault('RA_REF',[]).extend(GRa) 			#reference star RA
	matched_star_table.setdefault('Dec_REF',[]).extend(GDec)

	return


def distortion_plot_print_save(distortion_data,include,legendre_transformation,inverse_legendre_transformation,fp_iteration,ref_iteration,results):
	#this function takes a calculated distortion solution and generates various plots and outputs.
	#it needs to be passed info like the four-parameter iteration, reference iteration, year?
	#the main thing I want to save is the residuals.

	x = distortion_data['x_IM'].data
	y = distortion_data['y_IM'].data
	xref = distortion_data['x_REF'].data
	yref = distortion_data['y_REF'].data
	weights = distortion_data['Weight'].data
	star_id = distortion_data['Ref_id'].data
	xe = distortion_data['xe_IM'].data
	ye = distortion_data['ye_IM'].data
	xeref = distortion_data['xe_REF'].data
	yeref = distortion_data['ye_REF'].data
	frame = distortion_data['Frame_num'].data
	night_c = distortion_data['Night'].data
	ra = distortion_data['RA_IM'].data
	dec = distortion_data['Dec_IM'].data
	raref  = distortion_data['RA_REF'].data
	decref  = distortion_data['Dec_REF'].data	

	outliers = ~include

	xtf,ytf = legendre_transformation.evaluate(x,y)

	grid = np.arange(0,2048+1,64)
	x1,y1 = np.meshgrid(grid, grid)
	x2,y2 = legendre_transformation.evaluate(x1,y1)
	xgrid = x1.flatten()
	ygrid = y1.flatten()
	xgrid_tf = x2.flatten()
	ygrid_tf = y2.flatten()

	dr = np.sqrt((xgrid_tf-xgrid)**2 + (ygrid_tf-ygrid)**2)
	print("Mean distortion across field = ", np.mean(dr))

	if config['centred']:
		a, b = legendre_transformation.evaluate(frame_centre[0],frame_centre[1])
		x_central = a-frame_centre[0]
		y_central = b-frame_centre[1]
		xtf_c = xtf - x_central
		ytf_c = ytf - y_central
		xref_c = xref - x_central
		yref_c = yref - y_central
		xgrid_tf_c = xgrid_tf - x_central
		ygrid_tf_c = ygrid_tf - y_central
		dr = np.sqrt((xgrid_tf_c-xgrid)**2 + (ygrid_tf_c-ygrid)**2)
		print("Mean distortion after centering = ", np.mean(dr))


	residuals = np.hypot(xref-xtf,yref-ytf)     #distance between reference positions and corrected observed positions.

	# unexplained_x_error = 0.005
	# unexplained_y_error = 0.005
	unexplained_x_error = 0
	unexplained_y_error = 0
	unexplained_error = np.sqrt(unexplained_x_error**2 + unexplained_y_error**2)
	variances = 1/weights**2 + unexplained_error**2
	weights = 1/np.sqrt(variances)

	weighted_mean_residual = np.average(residuals[include], weights = weights[include])

	# print('min,median,max residuals:',np.min(residuals[include]),np.median(residuals[include]),np.max(residuals[include]))
	# print('min,median,max uncertainties:',np.min(1/weights[include]),np.median(1/weights[include]),np.max(1/weights[include]))
	# bad_a = np.nonzero(residuals > 22)
	# bad_o = np.nonzero(outliers) #manually finding outliers by their residuals
	# print('Residual_b > 22, not caught as outliers:', np.setdiff1d(bad_a, bad_o, assume_unique=True))

	#----------------------------------

	results.setdefault('median_residuals',[]).append(np.median(residuals[include]))
	results.setdefault('mean_residuals',[]).append(weighted_mean_residual)
	results.setdefault('min_residuals',[]).append(np.min(residuals[include]))
	results.setdefault('max_residuals',[]).append(np.max(residuals[include]))
	results.setdefault('num_residuals',[]).append(len(residuals[include]))

	with open(config['output_dir'] + 'iteration_residuals.txt'.format(), 'w') as temp:
		for r, resid in enumerate(results['median_residuals']):
			# temp.write(str(resid) + '\n')
			temp.write(f"Min:{results['min_residuals'][r]:.5f}  Median:{results['median_residuals'][r]:.5f}  Max:{results['max_residuals'][r]:.5f}  Num:{results['num_residuals'][r]:.5f}  Weighted mean:{results['mean_residuals'][r]:.5f}  \n") #Weighted mean squared:{results['mean_residuals_squared'][r]:.5f} Mean_a:{results['mean_residuals_r'][r]:.5f} Median_mas:{results['median_residuals_radec'][r]:.7f} 

	print('Ref iteration {} results:'.format(ref_iteration))
	print('Median residual = {}'.format(np.median(residuals[include])))
	print('Weighted mean residual = {}'.format(weighted_mean_residual))
	print('All ref_iteration median residuals = {}'.format(results['median_residuals']))
	print('All ref_iteration weighted mean residuals = {}'.format(results['mean_residuals']))


	#-----------plot difference vectors--------

	quiv_scale=200
	quiv_label_val = 10.0
	quiv_label = '{} pix'.format(quiv_label_val)
	margin = 100

	plt.figure(num=1,figsize=(6,6),clear=True)
	q = plt.quiver(xgrid,ygrid,(xgrid_tf-xgrid),(ygrid_tf-ygrid), color='black', scale=quiv_scale, angles='xy',width=0.003)
	plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
	plt.xlim(-margin,frame_size[0]+margin)
	plt.ylim(-margin,frame_size[1]+margin)
	# plt.axis('equal')
	# plt.set_aspect('equal','box')
	plt.title('Legendre polynomial model')
	plt.savefig('{}legendre_quiver_{}.jpg'.format(plot_dir,ref_iteration), bbox_inches='tight',dpi=200)

	plt.figure(num=1,figsize=(6,6),clear=True)
	q = plt.quiver(x[include],y[include],(xref[include]-x[include]),(yref[include]-y[include]), color='black', scale=quiv_scale, angles='xy',width=0.0005)
	q = plt.quiver(x[outliers],y[outliers],(xref[outliers]-x[outliers]),(yref[outliers]-y[outliers]), color='red', scale=quiv_scale, angles='xy',width=0.0005)
	plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
	plt.xlim(-margin,frame_size[0]+margin)
	plt.ylim(-margin,frame_size[1]+margin)
	# plt.axis('equal')
	# plt.set_aspect('equal','box')
	plt.title('Distortion of individual stars')
	plt.savefig('{}distortion_quiver_{}.jpg'.format(plot_dir,ref_iteration), bbox_inches='tight',dpi=200)


	if config['centred']:
		plt.figure(num=1,figsize=(6,6),clear=True)
		q = plt.quiver(xgrid,ygrid,(xgrid_tf_c-xgrid),(ygrid_tf_c-ygrid), color='black', scale=quiv_scale, angles='xy',width=0.003)
		plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
		plt.xlim(-margin,frame_size[0]+margin)
		plt.ylim(-margin,frame_size[1]+margin)
		# plt.axis('equal')
		# plt.set_aspect('equal','box')
		plt.title('Legendre polynomial model (centred)')
		plt.savefig('{}legendre_quiver_centred_{}.jpg'.format(plot_dir,ref_iteration), bbox_inches='tight',dpi=200)

		plt.figure(num=1,figsize=(6,6),clear=True)
		q = plt.quiver(x[include],y[include],(xref[include]-x_central-x[include]),(yref[include]-y_central-y[include]), color='black', scale=quiv_scale, angles='xy',width=0.0005)
		q = plt.quiver(x[outliers],y[outliers],(xref[outliers]-x_central-x[outliers]),(yref[outliers]-y_central-y[outliers]), color='red', scale=quiv_scale, angles='xy',width=0.0005)
		plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
		plt.xlim(-margin,frame_size[0]+margin)
		plt.ylim(-margin,frame_size[1]+margin)
		# plt.axis('equal')
		# plt.set_aspect('equal','box')
		plt.title('Distortion of individual stars (centred)')
		plt.savefig('{}distortion_quiver_centred_{}.jpg'.format(plot_dir,ref_iteration), bbox_inches='tight',dpi=200)


	#-----------------------------------------------------------------

	plt.figure(num=1,figsize=(6,6),clear=True)
	quiv_scale=5
	quiv_label_val = 0.5
	quiv_label = '{} pix'.format(quiv_label_val)
	q = plt.quiver(xtf[include],ytf[include],(xref[include]-xtf[include]),(yref[include]-ytf[include]),np.arctan2(yref[include]-ytf[include], xref[include]-xtf[include]),norm=Normalize(vmin=-math.pi,vmax=math.pi),  cmap='hsv', scale=quiv_scale, angles='xy',width=0.001)
	# q = plt.quiver(xtf[outliers],ytf[outliers],(xref[outliers]-xtf[outliers]),(yref[outliers]-ytf[outliers]), color='red', scale=quiv_scale, angles='xy',width=0.0005)
	plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
	plt.xlim(-margin,frame_size[0]+margin)
	plt.ylim(-margin,frame_size[1]+margin)
	# plt.axis('equal')
	# plt.set_aspect('equal','box')
	plt.title('Residual distortion'.format())
	# plt.savefig('{}/residual_b/residual_b_quiver_{}.pdf'.format(plot_dir,ref_iteration), bbox_inches='tight',dpi=200)
	plt.savefig('{}residual_quiver_{}.jpg'.format(plot_dir,ref_iteration), bbox_inches='tight',dpi=600)
	# plt.show()



def load_osiris_file(starfind_Dir,filename):
	#Loads a starlist file, converts it to a StarList object.
	print('Loading', filename)
	# with open (starfind_Dir + filename) as starfindFile:
			# starfindTable = starfindFile.read().splitlines()
	
	# osiris_starlist = starlists.StarList.from_lis_file(starfindTable, error=False)  #For .lis files


	osiris_table = Table.read(starfind_Dir + filename,format='ascii.basic',data_start=1)
	# osiris_table = Table(starfindTable)
	osiris_starlist = starlists.StarList.from_table(osiris_table)  				   #For _stars.txt files

	#from_lis_file requires these columns in order: name, m, t, x ,y , snr, corr, N_frames, flux
	# these are in the _stars.txt files, except for snr and corr.
	# modify my other script to generate appropriate files.
	#from_table requires name, x, y, and m.
	return osiris_starlist


def reference_section_plots(refTable_b,refTable_current):
	pass
	# plt.figure(num=1,figsize=(6,6),clear=True)
	# ptsize = 5
	# os.makedirs('{}combined_ref'.format(plot_dir),exist_ok=True)
	# # xbound = np.logical_and(-4.99 < refTable_b['x0'],refTable_b['x0'] < -4.77)
	# # ybound = np.logical_and(8.17 < refTable_b['y0'],refTable_b['y0'] < 8.4)
	# # s = np.where(np.logical_and(xbound,ybound))[0]
	# star_id = '86483'
	# s = np.where(refTable_b['name'] == star_id)[0]
	# plt.scatter(refTable_b['x0'][s],refTable_b['y0'][s],ptsize,'k',label='Mean')
	# plt.scatter(refTable_a['x'][s],refTable_a['y'][s],ptsize,label='Individual frames')
	# # plt.scatter(refTable_a['x'][s][0][0],refTable_a['y'][s][0][0],ptsize,label='Hubble?')
	
	# s = np.where(refTable_current[night]['name'].astype('str') == star_id)[0]
	# plt.scatter(refTable_current[night]['x0'][s],refTable_current[night]['y0'][s],ptsize,label='Ref')
	# plt.xlim([-4.932,-4.88])
	# plt.ylim([8.18,8.242])
	# plt.title('Averaging star position it:{}'.format(ref_iteration))
	# plt.legend(loc='upper left')
	# plt.xlabel('Relative RA (arcsec)')
	# plt.ylabel('Relative Dec (arcsec)')
	# plt.savefig('{}combined_ref/new_reference_stars_{}.jpg'.format(plot_dir,ref_iteration), bbox_inches='tight',dpi=200)
	# # plt.show()

	# print('{},{} matched stars in refTable_a'.format(len(idx1),len(idx2)))
	# #plot individual frame residuals
	# os.makedirs('{}residual_A/'.format(plot_dir,ref_iteration),exist_ok=True)
	# # for f in range(len(osiris_filenames)):
	# plt.figure(num=1,figsize=(6,6),clear=True)
	# # plt.clf()
	# quiv_scale=0.2
	# quiv_label_val = 0.01
	# quiv_label = '{} arcsec'.format(quiv_label_val)
	# # idx = np.where(frame[include] == f)[0]
	# # print(f'idx = {idx}')
	# # q = plt.quiver(xtf[include][idx],ytf[include][idx],(xref[include][idx]-xtf[include][idx]),(yref[include][idx]-ytf[include][idx]),np.arctan2(yref[include][idx]-ytf[include][idx], xref[include][idx]-xtf[include][idx]),norm=Normalize(vmin=-math.pi,vmax=math.pi),  cmap='hsv', scale=quiv_scale, angles='xy',width=0.005)
	# # q = plt.quiver(xtf[outliers],ytf[outliers],(xref[outliers]-xtf[outliers]),(yref[outliers]-ytf[outliers]), color='red', scale=quiv_scale, angles='xy',width=0.0005)
	# q = plt.quiver(refTable_b['x0'][idx1],refTable_b['y0'][idx1],(refTable_current[night]['x0'][idx2]-refTable_b['x0'][idx1]),(refTable_current[night]['y0'][idx2]-refTable_b['y0'][idx1]),np.arctan2(refTable_current[night]['y0'][idx2]-refTable_b['y0'][idx1], refTable_current[night]['x0'][idx2]-refTable_b['x0'][idx1]), norm=Normalize(vmin=-math.pi,vmax=math.pi), cmap='hsv', scale=quiv_scale, angles='xy',width=0.005)
	# plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
	# plt.xlim(-20,20)
	# plt.ylim(-20,20)
	# plt.xlabel('relative RA (arcsec)')
	# plt.ylabel('relative Dec (arcsec)')
	# # plt.axis('equal')
	# # plt.set_aspect('equal','box')
	# plt.title('Distortion residuals_r it:{}'.format(ref_iteration))
	# plt.savefig('{}residual_A/residual_quiver_{}.jpg'.format(plot_dir,ref_iteration), bbox_inches='tight',dpi=200)


def get_PA(filename):
	hdr = fits.getheader(filename)
	# pa = hdr['PA_IMAG']	
	pa = kai_instrument.get_instrument_angle(hdr)
	# return math.radians(hdr['PA_IMAG'])
	return pa

def calculate_PCU_guess(starlist_filename,night):
	if config['instrument'] == 'OSIRIS':
		flip_filename = starlist_filename[:-10] + '.fits'
		hdr = fits.getheader(config[night]['fits_dir']+flip_filename,ignore_missing_simple=True)
		
		if pcu_r_keyword not in hdr.keys():
			log = ascii.read(config[night]['log_file'],format='basic')
			raw_filename = flip_filename[:-10] + '.fits'
			mask = log['Filename']==raw_filename
			pcu_params = log[mask]
			hdr[pcu_x_keyword] = pcu_params['x'].data[0]
			hdr[pcu_y_keyword] = pcu_params['y'].data[0]
			hdr[pcu_z_keyword] = pcu_params['z'].data[0]
			hdr[pcu_r_keyword] = pcu_params['r'].data[0]

		pcu_x = float(hdr[pcu_x_keyword])   
		pcu_y = float(hdr[pcu_y_keyword])
		pcu_z = float(hdr[pcu_z_keyword])
		pcu_r = float(hdr[pcu_r_keyword])
		# pcu_r = 65.703

		#four parameter transformation.
		# x' = a0 + a1*x + a2*y
		# y' = b0 - a2*x + a1*y

		#in matrix notation:
		# x'  =  x  y 1 0  *  a1
		# y'  =  y -x 0 1     a2
		#					  a0
		#					  b0

		#from google:
		# a1 = S*cos(theta)
		# a2 = S*sin(theta)
		# so a2/a1 = tan(theta)
		# S = cos(theta)/a1

		# front_focus_scale = 1.385 #equals arcseconds per mm. 

		# #------Original, transforming pcu mm to pixels---------------------
		# S = 138.5  # 138.5 pixels per mm
		# theta_offset = 65.703   #this is some offset in degrees between the reported PCU rotation and the orientation we want.
		# theta = np.radians(pcu_r-theta_offset)  #check the sign of the rotation
		# a1 = S * math.cos(theta)
		# a2 = S * math.sin(theta)
		# pinhole_scale = 138.5 #pixles per mm on the pinhole mask.
		# pinhole_x_offset = 91.05  # mm
		# pinhole_y_offset = 183.95  # mm
		# pinhole_x_center = 1676 # pixels.  The location of the center of the pinohole mask in pixels when the PCU is at the offset position above.
		# # pinhole_y_center = 613 #pixels  On raw image
		# pinhole_y_center = 1422 #on flipped image
		# a0 = S*(pcu_x - pinhole_x_offset) + pinhole_x_center  #in pixels
		# b0 = -S*(pcu_y - pinhole_y_offset) + pinhole_y_center #flip applied to image, so using negative pix coordinates

		# #--------- pcu reference converted from mm to arcseconds, transformation from arcseconds to pixels--------------
		# S = 1/0.0099576 #pixels per arcsecond  This is used because the pinhole positions are converted to arcseconds before being passed to Flystar.
		# theta_offset = 65.703   #this is some offset in degrees between the reported PCU rotation and the orientation we want.
		# theta = np.radians(pcu_r-theta_offset)  #check the sign of the rotation
		# a1 = S * math.cos(theta)
		# a2 = S * math.sin(theta)
		# pinhole_scale = 138.5 #pixles per mm on the pinhole mask.
		# pinhole_x_offset = 91.05  # mm
		# pinhole_y_offset = 183.95  # mm
		# pinhole_x_center = 1676 #pixels.  The location of the center of the pinohole mask in pixels when the PCU is at the offset position above.
		# # pinhole_y_center = 613 #pixels  On raw image
		# pinhole_y_center = 1422 #on flipped image
		# a0 = S*(pcu_x - pinhole_x_offset) + pinhole_x_center  #in pixels
		# b0 = -S*(pcu_y - pinhole_y_offset) + pinhole_y_center #flip applied to image, so using negative pix coordinates


		#------forward transform: pixels to arcseconds------
		# S = 0.0099576 #arcseconds per pixel
		S = kai_instrument.get_plate_scale(hdr) #arcseconds per pixel
		pinhole_scale = 138.5 # pixels per mm
		front_focus_scale = S*pinhole_scale #arcsec per mm
		theta_offset = 65.703   #this is some offset in degrees between the reported PCU rotation and the orientation we want.
		theta = -np.radians(pcu_r-theta_offset)  #check the sign of the rotation
		a1 = S * math.cos(theta)
		a2 = S * math.sin(theta)
		pinhole_x_offset = 91.05  # mm
		pinhole_y_offset = 183.95  # mm
		# pinhole_x_center = 1676 #pixels.  The location of the center of the pinohole mask in pixels when the PCU is at the offset position above.
		# pinhole_y_center = 613 #pixels  On raw image
		# pinhole_y_center = 1422 #on flipped image
		pinhole_x_center = 1676 #pixels.  The location of the center of the pinohole mask in pixels when the PCU is at the offset position above.
		pinhole_y_center = 1422 #on flipped image
		a0 = -(S*pinhole_x_center + front_focus_scale*(pcu_x - pinhole_x_offset)) #in arcseconds
		b0 = -(S*pinhole_y_center + -front_focus_scale*(pcu_y - pinhole_y_offset)) #negative of inverse transformation (after scaling to arcsec)
		print('a0,b0',a0,b0)  #

		four_p = Empty()
		four_p.__class__ = transforms.four_paramNW   #I am trying to create an instance of the class without calling __init__, which requires two lists of stars to be passed in.
		four_p.px = [a0,a1,a2]
		four_p.py = [b0,-a2,a1]
		four_p.order = None
		four_p.mag_offset = 12 #need a way of calculating this.  Should be mean(m_ref-m_obs)
	else:
		print('Need to add PCU guess calculation for {}'.format(config['instrument']))
		sys.exit(0)
	return [four_p]

class Empty(object):
	#This exists so that we can create a transforms.four_paramNW object without passing it a list of star positions.
	pass


def mag_cut(starlist,min,max):
	idx = np.where((starlist['m'] >= min) & (starlist['m'] <= max))
	return starlist[idx]

def edge_cut(starlist,c,frame_size):
	x1 = c
	x2 = frame_size[0] - c 
	y1 = c 
	y2 = frame_size[1] - c 
	idx = np.where((starlist['x'] >= x1) & (starlist['x'] <= x2) & (starlist['y'] >= y1) & (starlist['y'] <= y2))
	# print('Cutting', len(starlist)-len(starlist[idx]), 'OSIRIS stars from the edge')
	return starlist[idx]

def trim_gaia(refTable,filename,PA):
	#select gaia points that are inside the frame
	with open (cleanDir + filename[:-12] + '.coo') as coofile:
		cooxy = np.genfromtxt(coofile)      #this won't work, I don't have coo files.
	# print('PA = ', PA)
	if abs(PA-0)<5:
		#S eclect a square, with x and y limits
		xmin=(-100 - cooxy[0]) * osiris_scale
		xmax=(frame_size[0]+100 - cooxy[0]) * osiris_scale
		ymin=(-100 - cooxy[1]) * osiris_scale
		ymax=(frame_size[1]+100 - cooxy[1]) * osiris_scale
		ind = (refTable['x0'] >= xmin) & (refTable['x0'] <= xmax) & (refTable['y0'] >= ymin) & (refTable['y0'] <= ymax)
	elif abs(PA-90)<5:
		xmax=-(-100 - cooxy[1]) * osiris_scale
		xmin=-(frame_size[0]+100 - cooxy[1]) * osiris_scale
		ymin=(-100 - cooxy[0]) * osiris_scale
		ymax=(frame_size[1]+100 - cooxy[0]) * osiris_scale
		ind = (refTable['x0'] >= xmin) & (refTable['x0'] <= xmax) & (refTable['y0'] >= ymin) & (refTable['y0'] <= ymax)
	else:
		# Select a circle centred on the frame
		x_centre = (frame_size[0]/2 - cooxy[0])*osiris_scale
		y_centre = (frame_size[1]/2 - cooxy[1])*osiris_scale
		# x_centre = 1024
		# y_centre = 1024
		# print(cooxy)
		# print(x_centre)
		# print(y_centre)
		# print(refTable['x0'])
		# print(refTable['x0']-x_centre)
		# x_centre = 0
		# y_centre = 0
		ind = np.hypot(refTable['x0']-x_centre, refTable['y0']-y_centre) <= np.hypot(1024,1024)*osiris_scale*2.0
		# print(ind)
		# quit()
	# print(np.hypot(1024,1024)*osiris_scale)
	# print(np.hypot(refTable['x0'], refTable['y0']))

	# print('Gaia bounds', xmin,xmax,ymin,ymax)
	# ind = (refTable['x0'] >= xmin) & (refTable['x0'] <= xmax) & (refTable['y0'] >= ymin) & (refTable['y0'] <= ymax)
	return refTable[ind]

	# wrong sometimes. 24, 27, Position angle

def plot_quiver(x,y,xref,yref,fitsfile,plot_directory,used,night,fp_iteration,ref_iteration):
	os.makedirs(plot_directory, exist_ok=True)		
	quiv_scale=100
	margin = 100
	# plt.close('all')
	plt.figure(num=4,figsize=(6,6),clear=True)
	# q = plt.quiver(x[~used],y[~used],(xref[~used]-x[~used]),(yref[~used]-y[~used]), color='red', scale=quiv_scale, angles='xy',width=0.004,label='Not UiT')
	q = plt.quiver(x[used],y[used],(xref[used]-x[used]),(yref[used]-y[used]), color='black', scale=quiv_scale, angles='xy',width=0.001)
	quiv_label = '10 pix'
	quiv_label_val = 10.0
	plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
	plt.xlim(-margin,frame_size[0]+margin)
	plt.ylim(-margin,frame_size[1]+margin)
	# plt.axis('equal')
	# plt.set_aspect('equal','box')
	plt.title('Individual star distortions, Observed -> Reference')
	# plt.savefig(plot_directory + fitsfile[:-5] + '_quiver.jpg', bbox_inches='tight')
	plt.savefig('{}{}_quiver_r{}_f{}.jpg'.format(plot_directory,fitsfile[:-5],ref_iteration,fp_iteration), bbox_inches='tight',dpi=200)

	# plt.show()
	# plt.close()

def plot_scatter(x,y,fitsfile,plot_directory,night,fp_iteration,ref_iteration):
	os.makedirs(plot_directory, exist_ok=True)		
	# fitsfile = cleanDir + filename[:-12] + '.fits'
	img = fits.getdata(config[night]['fits_dir']+fitsfile)
	img[np.where(img<1)] = 1
	# plt.close('all')
	plt.figure(num=4,figsize=(6,6),clear=True)
	margin = 100
	vmin = 100
	vmax = np.max(img)*0.75
	norm = LogNorm(vmin, vmax)
	plt.imshow(img, cmap='Greys_r', norm=norm, origin = "lower", )
	ptsize = 10
	sc = plt.scatter(x,y,s=ptsize,c='g',alpha=1.0,marker='.',label='Observed')
	plt.xlim(-margin,frame_size[0]+margin)
	plt.ylim(-margin,frame_size[1]+margin)
	plt.title('Observed stars')
	plt.legend(loc='upper right')
	# plt.savefig(plot_directory + fitsfile[:-5] + '_scatter.pdf', bbox_inches='tight')
	plt.savefig('{}{}_scatter_r{}_f{}.jpg'.format(plot_directory,fitsfile[:-5],ref_iteration,fp_iteration),bbox_inches='tight',dpi=200)

def plot_matched(x,y,xref,yref,fitsfile,plot_directory,used,night,fp_iteration,ref_iteration):
	os.makedirs(plot_directory, exist_ok=True)		
	# fitsfile = cleanDir + filename[:-12] + '.fits'
	img = fits.getdata(config[night]['fits_dir']+fitsfile)
	img[np.where(img<1)] = 1
	# plt.close('all')
	plt.figure(num=4,figsize=(6,6),clear=True)
	margin = 100
	vmin = 100
	vmax = np.max(img)*0.75
	norm = LogNorm(vmin, vmax)
	plt.imshow(img, cmap='Greys_r', norm=norm, origin = "lower", )
	ptsize = 10
	plt.scatter(xref,yref,s=ptsize,c='r',alpha=0.6,marker='+',label='Reference')
	# sc = plt.scatter(x[~used],y[~used],s=ptsize,c='b',marker='.',label='Observed not UiT')
	sc = plt.scatter(x[used],y[used],s=ptsize,c='g',alpha=0.6,marker='.',label='Observed')
	plt.xlim(-margin,frame_size[0]+margin)
	plt.ylim(-margin,frame_size[1]+margin)
	plt.title('Matched stars')
	plt.legend(loc='upper right')
	# plt.savefig(plot_directory + fitsfile[:-5] + '_matched.pdf', bbox_inches='tight')
	plt.savefig('{}{}_matched_r{}_f{}.jpg'.format(plot_directory,fitsfile[:-5],ref_iteration,fp_iteration), bbox_inches='tight',dpi=200)


def distortion_from_coefficients(filename):
	#reads a file with a list of coefficients, generates a transformation object
	if os.path.exists(filename):
		coefficients = ascii.read(filename,format='fixed_width')
		px = coefficients['px_val']
		py = coefficients['py_val']
		order = int(np.sqrt(len(px))-1)
		x_domain = [0,frame_size[0]] #I am not sure if this is the correct domain. Could be [0,1]? It doesn't actually get used anywhere, so fine for now.
		y_domain = [0,frame_size[1]]
		tform = transforms.LegTransform(order, px, py, x_domain, y_domain, astropy_order=True)
	else:
		print('Previous distortion model not found: {} '.format(filename))
		tform = None
	return tform



def stack_frames():
	pass



def load_reference_list():
	#remember, there can be different lists in one run
	pass


def load_observed_lists():
	pass



#--------------------------------
main()
sys.exit(0)
#--------------------------------







