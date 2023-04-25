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

command looks like: python bdc_main.py config1.yml


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
matplotlib.use('agg') #stops memory leak?
from matplotlib.colors import LogNorm, Normalize
import numpy as np
import math
import datetime

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config", required = True, help = "Configuration file")
parser.add_argument("-o", "--output", required = True, help = "Output folder")
args = parser.parse_args()
print("Config file is: {}".format(args.config))
print("Output folder is : {}".format(args.output))

with open(args.config, "r") as ymlfile:
	config = yaml.safe_load(ymlfile)

	# for section in config:
	# 	print(section)
	# print('----')
	# # print(cfg)
	# print(config["reference"])
	# print(config["n1"]['bad_files'])
	# print(config["n1"]['bad_files'][0])
	# for night in config['nights']:
	# 	print(config[night]['target'])


def main():
	if config['instrument'] == 'OSIRIS'
		osiris = instruments.OSIRIS()
	else:
		print('Can\'t handle instrument {}, only OSIRIS'.format(config['instrument']))	
	
	refTable_current = {}	#holds most recent refTable for each night #this is the reference table that keeps getting updated by the reference section.
	if config['generate_reference_frame']:
		for ref_iteration in range(config['ref_iterations']):
			print('\n \n Ref Iteration {} \n'.format(ref_iteration))
			distortion_model = distortion_section(refTable_current)
			new_refTable = reference_section(refTable_current,distortion_model)
			refTable_current = new_refTable
	else:
		distortion_model = distortion_section(refTable_current)
		#outputs are printed to files.


	print('Completed ref_iterations')
	#here is outside the iteration loop
	print('Median residuals_b for each ref_iteration')
	print(median_residuals_b)
	print('Max residuals_b for each ref_iteration')
	print(max_residuals_b)
	print('Min residuals_b for each ref_iteration')
	print(min_residuals_b)
	print('Number of residuals_b for each ref_iteration')
	print(num_residuals_b)
	print('Weighted mean residuals_b for each ref_iteration')
	print(mean_residuals_b)
	print('Weighted mean residuals_b squared for each ref_iteration')
	print(mean_residuals_b_squared)

	return distortion_model

def distortion_section(refTable_current):
	#this section takes a reference table and some OSIRIS observations, matches them, and calculates a distortion model using Legendre polynomials.

	# current_distortion_correction = initial_distortion_correction
	current_distortion_correction = None
	tab1_initial = {} #we only want to use the match table from the first four-parameter iteration
	for night in obs_nights:
		tab1_initial[night]=[]

	refTable_current_filename = '{}refTable_current_{}_{}'.format(resultDir,fp_iteration) #this should be somewhere else?

	# do I need the empty lists here in order to append to them?
	#if I am loading dist files then tab1_intial can break if a set of fp_iterations was incomplete - it won't run the first iteration to generate tab1_initial. Need to save as a file?

	for fp_iteration in fp_iterations:

		matched_star_table_name = '{}dist_measures_{}_{}_{}.txt'.format(resultDir,solution_year,ref_iteration,fp_iteration)
		matched_star_table = {}   #matched star lists from each image are all appended to this dictionary. Should be an astropy table? Lots of appending lists. Convert later?
		plate_scale_and_rotations = {}


		if config['generate_new_fp']:  		#if config['generate_new_fp']==True, we generate a new table of matched stars. If False, we load the previouly saved one.
			obs_nights = [night for night in config]
			#-------------------------Night loop starts here-------------------------------
			for night in obs_nights:
				hubbleData = fetch_hubble(config[night]['reference_file'])
				osiris_filenames = get_image_filenames(config[night]['starlist_dir'],config[night]['slice'])
				if night not in refTable_current.keys():  #if we have not generated a new reference frame, load the starting one.
					if config['use_flystar_velocity']:   #if True, use flystar to project positions. If False, manually project positions.
						refTable_H = prepare_ref_for_flystar(hubbleData,config[night]['target_RA'],config[night]['target_DEC'],config[night]['target'],config[night]['ref_instrument'])
					else:
						starlist0 = load_osiris_file(config[night]['starlist_dir'] ,osiris_filenames[0]) #load 1 file to get the epoch.
						hubbleData_p = project_pos(hubbleData,starlist0,'hubble_'+config[night]['target'])
						refTable_H = prepare_ref_for_flystar(hubbleData_p,config[night]['target_RA'],config[night]['target_DEC'],config[night]['target'],config[night]['ref_instrument'])
					refTable_current[night] = refTable_H.filled() #is this necessary?
					with open(refTable_current_filename, 'wb') as temp:
						pickle.dump(refTable_current, temp)

				for i, filename in enumerate(osiris_filenames):
					if filename in config[night]['bad_files']:
						print('{} {} flagged as bad, skipping'.format(i,filename))
						continue

					starlist = load_osiris_file(config[night]['starlist_dir'] ,filename)
					fitsfile = config[night]['cleanDir'] + filename[:-12] + '.fits'
					PA = get_PA(fitsfile)
					# starlist = brightest_n(starlist,170)
					starlist = mag_cut(starlist,0,config[night]['minmag'])
					if not filename in config[night]['dont_trim']:
						starlist = edge_cut(starlist,5)
					if len(starlist) == 0:					
						print(i,filename, '0 stars remaining after edge cut, skipping image')     #this needs to be outside the function.
						errorlist.append(filename[:-4])
						continue
					refTable_t = trim_gaia(refTable_current[night],filename,PA)    #this will need to be changed for the pinhole mask.
					refTable_tf = refTable_t.filled()		#unmask, required for quiver plots. Can this be removed?
					if config['ref_instrument'] == 'PCU':
						refTable_d = refTable_tf[:]   #The PCU does not need to have DAR applied
					else:
						refTable_d = dar.applyDAR(fitsfile, refTable_tf, plot=False, instrument=osiris, plotdir=plotDir_n + 'dar_a/'+ str(ref_iteration) + '/')

					try:

						transformation_guess_file = '{}tform_{}_{}_{}.p'.format(tformDir,ref_iteration,fp_iteration,filename) #this needs to be reworked for the PCU

						if config['ref_instrument'] == 'PCU':
							transformation_guess = calculate_PCU_guess(filename)
						else:
							if os.path.exists(transformation_guess_file):
								print('Loading transform guess')
								with open(transformation_guess_file, 'rb') as trans_file:
									transformation_guess = pickle.load(trans_file) 	
							else:
								# transformation_guess = last_good_transform
								print('Not loading transform guess')
								transformation_guess = None


						starlist_corrected = starlist[:]
						# if current_distortion_correction is not None:    #either this check, or the fp_iteration check.
						if fp_iteration > 0:
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
							    init_guess_mode='miracle', verbose=0)
							msc.fit()
							# tab1 = msc.ref_table
							tform = msc.trans_list
							tform_inv = msc.trans_list_inverse
							#I only want to save the transformation from the subsequent iterations, not tab1.
							if fp_iteration == 0:
								tab1_initial[night].append(msc.ref_table)
								#save tab1_initial to a file. Once it has all nights in it?
							tab1 = tab1_initial[night][i]
							#
							j = 0
						else:
							j=i  # i=osiris index, j=gaia index
							print('Only set up for single_fit mode')

						ref_idx = np.where(tab1['ref_orig'] == True)[j]
						if len(ref_idx) <= 10:
							print(i, filename, 'Only', len(ref_idx),'matched, skipping')
							errorlist.append(filename[:-4])
							return
						with open(transformation_guess_file, 'wb') as temp:
							pickle.dump(tform, temp)

						print(i, filename, len(refTable_d), 'Reference stars,', len(starlist_corrected), 'OSIRIS stars,', len(ref_idx), 'matches')
						
						px = tform_inv[j].px
						py = tform_inv[j].py
						theta = math.atan2(px[2],px[1])
						scale = math.cos(theta) / px[1]
						#this may be a mistake, google says it's px[1]/cos(theta). Maybe that's why the inverse transformation works?
						plate_scale_and_rotations.setdefault('Filename',[]).append(filename)
						plate_scale_and_rotations.setdefault('Scale',[]).append(scale)
						plate_scale_and_rotations.setdefault('Rotation',[]).append(math.degrees(theta))
						plate_scale_and_rotations.setdefault('PA',[]).append(PA)

						append_to_matched_star_table(matched_star_table,tab1,tform,tform_inv)

						# plot_quiver(tab1['x0'][ref_idx],tab1['y0'][ref_idx],tab1['x'][ref_idx,0],tab1['y'][ref_idx,0],filename[:-4])

					except AssertionError as err:
						print(filename[:-4], 'error:')
						print(err)
						errorlist.append(filename[:-4])
						continue

					except ValueError as err:
						print(filename[:-4], 'error:')
						print(err)
						errorlist.append(filename[:-4])
						continue

					successlist.append(filename[:-4])

				#--------------------------------------------Night loop finishes here----------------------------------------------

				print('--------------------')
				print('Succeeded for', len(successlist), 'files.')
				print('Failed for', len(errorlist), 'files:')
				print(errorlist)

				plate_scale_and_rotations['Difference'] = plate_scale_and_rotations['Rotation'] - plate_scale_and_rotations['PA']
				transform_table = Table(plate_scale_and_rotations)
				ascii.write(transform_table,'{}t_params_{}_{}.txt'.format(resultDir,ref_iteration,night),format='fixed_width', overwrite=True)
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

		outliers, bad_names, m_distances = find_outliers(distortion_data['x_IM'],distortion_data['y_IM'],distortion_data['x_REF'],distortion_data['y_REF'],distortion_data['REF_ID'])
		include = ~outliers 
		print('Mahalanobis distance. Mean: {}, STD: {}'.format(np.mean(m_distances),np.std(m_distances)))

		print('Fitting Legendre Polynomial')
		legendre_transformation = transforms.LegTransform.derive_transform(distortion_data['x_IM'][include], 
																			distortion_data['y_IM'][include], 
																			distortion_data['x_REF'][include], 
																			distortion_data['y_REF'][include], 
																			config['legendre_order'], m=None, mref=None,init_gx=None, init_gy=None, weights=None, mag_trans=True
																			) # Defines a bivariate legendre tranformation from x,y -> xref,yref using Legnedre polynomials as the basis.

		current_distortion_correction = legendre_transformation  #update current_distortion_correction with latest Legendre polynomial

		distortion_plot_print_save(distortion_data,current_distortion_correction,fp_iteration)

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
		ascii.write(output_table,'{}distortion_coefficients_{}_{}_{}.txt'.format(resultDir,solution_year,ref_iteration,fp_iteration),format='fixed_width', overwrite=True)


	return current_distortion_correction #returns the final distortion correction




def reference_section(refTable_current,distortion_model):
	#this function applies a distortion solution to the OSIRIS lists, and combines them into a single ref list to return (per night)
		if len(refTable_current) == 0:
			with open(refTable_current_filename, 'rb') as temp:
				refTable_current = pickle.load(temp)   #Currently I do need this, to load the reftable saved from ... the last time.

		#------------------------Night loop 1 starts here---------------------------

		obs_nights = [night for night in config]
		print(obs_nights)
		for night in obs_nights:

			plotDir_n = plotDir + night + '/'   #Need to sort out which plots I actually need.

			osiris_filenames = get_image_filenames(config[night]['starlist_dir'],config[night]['slice'])
			print(len(osiris_filenames), 'OSIRIS images')

			#perhaps make this next section possibly generate a new reference frame, then always load from file. For consistency.
			#make a function to choose the ref file? It's just checking if there's a previous one, can probably do one step.

			# combined_ref_filename = resultDir + 'combined_ref_table_' + str(ref_iteration) + '.txt'
			combined_ref_filename = '{}combined_ref_table_{}_{}.txt'.format(resultDir,night,ref_iteration)

			if create_combined_reflist or (not os.path.exists(combined_ref_filename)):
				list_of_starlists = []
				print('Generating new combined refTable')

				for i, filename in enumerate(osiris_filenames):
					if filename in config[night]['bad_files']:
						print('{} {} flagged as bad, skipping'.format(i,filename))
						continue
						# errorlist.append(filename[:-4])
					print('{} {} applying distortion correction'.format(i,filename))
					
					fitsfile = config[night]['cleanDir'] + filename[:-12] + '.fits'
					PA = get_PA(fitsfile)
					starlist = load_osiris_file(config[night]['config[night]['starlist_dir']'] ,filename)
					# starlist = brightest_n(starlist,170)
					starlist = mag_cut(starlist,0,minmag)
					if not filename in config[night]['dont_trim']:
						starlist = edge_cut(starlist,5)
					if len(starlist) == 0:
						print(i,filename, '0 stars remaining after edge cut, skipping image')
						errorlist.append(filename[:-4])
						continue

					#load_and_prepare_data does something similar, but also trims the reference. Maybe split that up and use the function.
					#or just don't use the function, it's only a few lines.

					xt, yt = correction_list[ref_iteration].evaluate(starlist['x'],starlist['y'])    #apply distortion correction

					if config[night]['centred']:
						xc, yc = correction_list[ref_iteration].evaluate(1024,1024)
						xc -= 1024
						yc -= 1024						
					else:
						xc = 0
						yc = 0
					starlist['x'] = xt - xc
					starlist['y'] = yt - yc

					# plt.figure(num=4,figsize=(6,6),clear=True)			
					if config['ref_instrument'] != 'PCU':
						starlist = dar.removeDAR(fitsfile,starlist, plot=False, instrument=osiris, plotdir=plotDir_n + 'dar_r/'+ str(ref_iteration) + '/')

					list_of_starlists.append(starlist)

					# successlist.append(filename[:-4])

				print('--------------------')
				print('Completed applying distortion corrections')

		
				reference_transformation_guess_file = '{}combining_ref/tform_{}_{}.p'.format(tformDir,ref_iteration,night)

				# this should be a list of all the tranformation guesses from the distortion section.
				single_transformation_guess_list = ['/u/mfreeman/work/d/transform_files/hubble/tform_{}.p'.format(i) for i in osiris_filenames]

				use_single_trans_guesses = False
				if os.path.exists(reference_transformation_guess_file):
					print('Loading transform from {}'.format(reference_transformation_guess_file))
					with open(reference_transformation_guess_file, 'rb') as trans_file:
						trans_list_temp = pickle.load(trans_file)

					if len(trans_list_temp) == len(list_of_starlists):
						print('Transform guess has correct number of lists {}'.format(len(trans_list_temp)))
						reference_transformation_guess = trans_list_temp
					else:
						use_single_trans_guesses = True
				else:
					use_single_trans_guesses = True

				if use_single_trans_guesses:
					print('Loading old transform')
					if os.path.exists(tform_file_ref_previous[0]):
						reference_transformation_guess = []
						for i, tform_filename in enumerate(tform_file_ref_previous):
							with open(tform_filename, 'rb') as trans_file:
								trans_list.extend(pickle.load(trans_file))					#initial guess for the transformation
																							#each file is a list with one element, so using .extend
					else:
						# trans_list = last_good_transform
						print('No trans_list found')
						reference_transformation_guess = None


				#Do I need to be able to pass these parameters in, or should they stay fixed? My testing found False/False to work best.
				set_use_ref_new = False
				set_update_ref_orig = False

				print('Creating combined reference frame with MosaicToRef')
				msc = align.MosaicToRef(refTable_current[night],
					list_of_starlists, iters=3,
					dr_tol=[0.5,0.5,0.5], dm_tol=[2,2,2],
					outlier_tol=[None,None,None],
					trans_input=reference_transformation_guess,
					trans_class=transforms.four_paramNW,
					trans_args=[{'order': 1}, {'order': 1}, {'order': 1}],
					use_vel=config['use_flystar_velocity'],
					mag_trans=True,
					mag_lim=[6,13], #[6,16],
					weights=None,
					use_ref_new = set_use_ref_new,
					update_ref_orig = set_update_ref_orig,
					calc_trans_inverse=True,    
					init_guess_mode='miracle', verbose=0)


				msc.fit()
				# tab2 = msc.ref_table
				with open(reference_transformation_guess_file, 'wb') as temp:
					pickle.dump(msc.reference_transformation_guess, temp)
				print('Completed reference table matching')
				# refTable_a = msc.ref_table
				with open(combined_ref_filename, 'wb') as temp:
					pickle.dump(msc.ref_table, temp)

			else:
				print('a: Loading previous combined refTable')

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

			#can I just copy the whole table and delete the 'use_in_trans' column? The current way works at least.

			if trim_not_used_in_trans:
				print('Trimming refTable_b to only include stars used in transformation')
				refTable_b = refTable_b[refTable_a['use_in_trans']]			#using refTable_a due to the problem above with refTable_b

			temp, idx1, idx2 = np.intersect1d(refTable_b['name'].astype('str'),refTable_current[night]['name'],return_indices=True)

			reference_section_plots(refTable_b)

			residuals_a = np.hypot(refTable_current[night]['x0'][idx2]-refTable_b['x0'][idx1],refTable_current[night]['y0'][idx2]-refTable_b['y0'][idx1])
			mean_residuals_a.append(np.mean(residuals_a))  #I don't want the hypot, keep them separate? Should be the mean anyway.

			refTable_current[night] = Table(refTable_b,copy=True)

	return refTable_current








def find_outliers(x,y,xref,yref,star_id,verbose=False):
	print('Finding outliers')
	#scan grid over field. Flag any points as outliers.  Service(2016) did 205x205 pixel bins. removing 3 sigma outliers.
	dx = x - xref
	dy = y - yref
	# bins = np.array([0,512,1024,1536,2048])   #looks like 512x512 bins are the best.
	bins = np.arange(0,2048+1,256)

	xbin = np.digitize(x,bins)
	ybin = np.digitize(y,bins)

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
			for i, xb in enumerate(bins):
				for j, yb in enumerate(bins):
					
					if i == 0 or j == 0:
						continue

					include = ~outflag
					idx = np.nonzero((xbin[include] == i) & (ybin[include] == j))[0]
					
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
		for i, xb in enumerate(bins):
			for j, yb in enumerate(bins):
				
				if i == 0 or j == 0:
					continue
				
				include = ~outflag
					# print(idx)
				idx = np.nonzero((xbin == i) & (ybin == j))[0]


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





def fetch_hubble(filename):				#loads hubble data from hst1pass
	# filelist = os.listdir(nightDir)
	# hst1pass_file = fnmatch.filter(filelist,'i*.xymrduvqpk')[0]
	print('Loading', filename)
	# hubble_table = ascii.read(nightDir + hst1pass_file,format='commented_header',guess=False,delimiter='\s',comment='#',header_start=166,data_start=168,) #lines start at 0
	# hubble_table = ascii.read(nightDir + hst1pass_file,format='basic',names=['x','y','m','r','d','u','v','q','p','k'])
	if instrument == 'Hubble'
		column_names = ['r', 'x_0', 'y_0', 'pmx', 'pmy', '1sig_pmx', '1sig_pmy', 'x_M', 'y_M', 'Delta_x', 'Delta_y', 
						'1sig_pmx_mod', '1sig_pmy_mod', 'm_F606W', 'm_F814W', 'rms_F606W', 'rms_F814W', 'qfit_F606W', 'qfit_F814W', 
						'red_chi2_pmx', 'red_chi2_pmy', '1sig_intrcp_pmx', '1sig_intrcp_pmy', 'time_bsln', '1sig_intrcp_pmx_mod', '1sig_intrcp_pmy_mod', 
						'pm_ref', 'N_found', 'N_used', 'ID', 'delta_pmx', 'delta_pmy']
		hubble_table = ascii.read(filename,format='basic',names=column_names)
	
	if instrument == 'PCU'
		hubble_table = ascii.read(filename,format='fixed_width')

	return hubble_table

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
		front_focus_scale = 1.38 #milliarcseconds per micron, which equals arcseconds per mm.  
		#Measurements with pcu give 138.5 pixels per mm -> 1385 mas per mm -> 1.385 arcsec per mm.


		ref_new = Table([ref['id']], names=['name'], masked=False)

		ref_new['x0'] = (ref['x'] * front_focus_scale) *-1            #*-1 to get East positive
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


def get_image_filenames(directory,file_slice):
	filelist = os.listdir(directory)
	# starfindFiles = fnmatch.filter(filelist,'ci*.lis')
	image_files = fnmatch.filter(filelist,'i*_stars.txt')  #these may need to be converted to .lis files?
	image_files.sort()
	sl = slice(file_slice[0],file_slice[1])
	image_files = image_files[sl]	
	return image_files

def project_pos(refData_o, starlist,instrument):
	#Assume that all osiris images have the same epoch.
	#GAIA['ra'] and GAIA['dec'] are in degrees, GAIA['pmra'] and GAIA['pmdec'] are in mas/year
	#hubble['x_0'] is in arcseconds, hubble['pmx'] is in mas per year
	refData = refData_o[:]   #make a copy of the original array
	osirisT = starlist['t'][0]
	if instrument == 'gaia':
		dt = osirisT - refData['ref_epoch']
		print('Projecting Gaia {:.3f} years'.format(dt[0]))
		idra=~refData['pmra'].mask    #select not masked
		iddec=~refData['pmdec'].mask 	
		refData['ra'][idra] = refData['ra'][idra] + refData['pmra'][idra]/(3600*1000) * dt[idra]        #mas to degrees
		refData['dec'][iddec] = refData['dec'][iddec] + refData['pmdec'][iddec]/(3600*1000) * dt[iddec]
		#need to include error, like hubble below.
		refData['ref_epoch'] = osirisT  #??

	elif instrument[0:6] =='hubble':
		if instrument[-3:] == 'm15':
			#m15 reference time is J2006.334577
			hubbleT = 2006.334577 								#hard coded m15
		elif instrument[-3:] == 'm92':
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

	else:	 
		print(instrument, 'in project_pos should be hubblem## or gaia')
	return refData	


def load_transformation_guess():
	if reference_instrument == 'Hubble':
		tform_file_1 = '{}tform_{}_{}_{}.p'.format(tformDir,ref_iteration,fp_iteration,filename)	#don't need to enter the night, because the filename includes the date and so is unique
		tform_file_2 = '{}tform_{}_{}.p'.format(transform_files_location,ref_iteration,filename)
		tform_file_3 = '/u/mfreeman/work/d/transform_files/hubble/tform_{}.p'.format(filename)

	elif reference_instrument == 'GAIA':
		#out of date
		tform_file_1 = '{}tform_{}_{}.p'.format(tformDir,ref_iteration,filename)
		tform_file_2 = '{}tform_{}_{}.p'.format(transform_files_location,fitmode,ref_iteration,filename)
		tform_file_3 = '/u/mfreeman/work/d/transform_files/gaia/tform_{}.p'.format(filename)

	if os.path.exists(tform_file_1):
		print('Loading correct transform guess')
		with open(tform_file_1, 'rb') as trans_file:
			trans_list = pickle.load(trans_file) 					#initial guess for the transformation
	elif os.path.exists(tform_file_2):
		print('Loading default transform guess')
		with open(tform_file_2, 'rb') as trans_file:
			trans_list = pickle.load(trans_file)
	elif os.path.exists(tform_file_3):
		print('Loading old transform guess')
		with open(tform_file_3, 'rb') as trans_file:
			trans_list = pickle.load(trans_file)
	else:
		# trans_list = last_good_transform
		print('Not loading transform guess')
		trans_list = None	
	return trans_list


def append_to_matched_star_table(matched_star_table,tab1,ref_idx,tform,tform_inv);
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

	#don't need these plots?
	plot_image(Gx,Gy,Ox,Oy,fitsfile,plotDir_n + 'img/'+ str(ref_iteration) + '/',tab1['use_in_trans'][ref_idx])
	plot_image_dots(Gx,Gy,Ox,Oy,fitsfile,plotDir_n + 'img_d/'+ str(ref_iteration) + '/',tab1['use_in_trans'][ref_idx])
	plot_quiver(Ox,Oy,Gx,Gy,filename[:-4],plotDir_n + 'quiver/'+ str(ref_iteration) + '/',tab1['use_in_trans'][ref_idx])
	if show_plots:
		plt.show()

	return


def distortion_plot_print_save(distortion_data,legendre_transformation,fp_iteration):
		#this function takes a calculated distortion solution and generates various plots and outputs.
		#it needs to be passed info like the four-parameter iteration, reference iteration, year?
		#the main thing I want to save is the residuals.


		x = distortion_data['x_IM'].data
		y = distortion_data['y_IM'].data
		xref = distortion_data['x_REF'].data
		yref = distortion_data['y_REF'].data
		weights = distortion_data['Weight'].data
		star_id = distortion_data['REF_ID'].data
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


		xts,yts = legendre_transformation.evaluate(x,y)
		# unexplained_x_error = 0.005
		# unexplained_y_error = 0.005
		unexplained_x_error = 0
		unexplained_y_error = 0
		if centred:
			a,b = legendre_transformation.evaluate(1024,1024)
			x_central = a-1024
			y_central = b-1024
			xts = xts - x_central
			yts = yts - y_central
			xref = xref - x_central
			yref = yref - y_central

		xts_ra = []
		yts_dec = []

		#-------------------plotting quivers in RA Dec ------------------- 


		for night in obs_nights:
			# plt.figure(num=1,figsize=(6,6),clear=True)
			# plt.clf()
			# plt.figure(num=2,figsize=(6,6),clear=True)
			# plt.clf()
			quiv_scale= 0.05
			quiv_label_val = 0.001
			quiv_label = '{} mas'.format(quiv_label_val*1000)

			os.makedirs('{}sky_resid_b_individual/{}/{}/{}'.format(plotDir,ref_iteration,night,fp_iteration),exist_ok=True)

			for f, filename_1 in enumerate(osiris_filenames_dict[night]):
				# tform_file_4p = './transform_files/hubble/tform_{}_{}.p'.format(ref_iteration,filename_1)
				tform_file_4p = '{}tform_{}_{}_{}.p'.format(tformDir,ref_iteration,fp_iteration,filename_1)
				with open(tform_file_4p, 'rb') as trans_file:
					transform_4p = pickle.load(trans_file) 	

				# idx = np.where(frame[include] == f)[0]
				idx = np.where((frame[include] == f) & (night_c[include] == night))[0]
				# RA_ref, Dec_ref = transform_4p[0].evaluate(xref[include][idx],yref[include][idx])	#transform pixel coordinates to RA/Dec
				# RA_xts, Dec_yts = transform_4p[0].evaluate(xts[include][idx],yts[include][idx])  #distortion corrected, converted to RA/Dec

				RA_ref = raref[include][idx]	#instead of applying the inverse transformation then the transformation, just use the orignal ref coords. Not Centred
				Dec_ref = decref[include][idx]
				RA_xts, Dec_yts = transform_4p[0].evaluate(xts[include][idx]+ x_central,yts[include][idx]+y_central) #not centred
				xts_ra.extend(RA_xts)
				yts_dec.extend(Dec_yts)

				# print(f'idx = {idx}')
				angle_colour = np.arctan2(Dec_ref-Dec_yts, RA_ref-RA_xts)
				# norm = Normalize(vmin=0,vmax=2*math.pi)
				# norm.autoscale(angle_colour)
				# colourmap = 'hsv'
				# q = plt.quiver(RA_xts,Dec_yts,(RA_ref-RA_xts),(Dec_ref-Dec_yts),angle_colour, norm=Normalize(vmin=-math.pi,vmax=math.pi), cmap='hsv', scale=quiv_scale, angles='xy',width=0.002)	#from corrected to ref
				plt.figure(num=2,figsize=(6,6),clear=False)
				q = plt.quiver(RA_ref,Dec_ref,(-RA_ref+RA_xts),(-Dec_ref+Dec_yts),angle_colour, norm=Normalize(vmin=-math.pi,vmax=math.pi), cmap='hsv', scale=quiv_scale, angles='xy',width=0.002)  #from ref to corrected
				
				plt.figure(num=1,figsize=(6,6),clear=True)
				# plt.clf()
				# q = plt.quiver(RA_xts,Dec_yts,(RA_ref-RA_xts),(Dec_ref-Dec_yts),angle_colour, norm=Normalize(vmin=-math.pi,vmax=math.pi), cmap='hsv', scale=quiv_scale, angles='xy',width=0.002)
				q = plt.quiver(RA_ref,Dec_ref,(-RA_ref+RA_xts),(-Dec_ref+Dec_yts),angle_colour, norm=Normalize(vmin=-math.pi,vmax=math.pi), cmap='hsv', scale=quiv_scale, angles='xy',width=0.002)
				plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
				plt.xlim(-20,20)
				plt.ylim(-20,20)
				plt.xlabel('RA (arcsec)')
				plt.ylabel('Dec (arcsec)')
				plt.title('Sky Distortion residuals_b {} {}'.format(ref_iteration,f))
				plt.savefig('{}sky_resid_b_individual/{}/{}/{}/residual_b_quiver_{}_{}_{}_{}_{}.pdf'.format(plotDir,ref_iteration,night,fp_iteration,solution_year,ref_iteration,night,fp_iteration,f), bbox_inches='tight',dpi=200)
				# q = plt.quiver(xts[outliers],yts[outliers],(xref[outliers]-xts[outliers]),(yref[outliers]-yts[outliers]), color='red', scale=quiv_scale, angles='xy',width=0.0005)

			plt.figure(num=2,figsize=(6,6),clear=False)
			plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
			plt.xlim(-20,20)
			plt.ylim(-20,20)
			plt.xlabel('RA (arcsec)')
			plt.ylabel('Dec (arcsec)')
			plt.title('Sky Distortion residuals_b {}'.format(ref_iteration))
			os.makedirs('{}sky_resid_b/{}/{}/'.format(plotDir,ref_iteration,night),exist_ok=True)
			plt.savefig('{}sky_resid_b/{}/{}/residual_b_quiver_{}_{}_{}_{}.pdf'.format(plotDir,ref_iteration,night,solution_year,ref_iteration,night,fp_iteration), bbox_inches='tight',dpi=200)
		# plt.close('all')

		#-------------------------------------------------------------
		#------------------Plot quivers in pixels--------------------
		os.makedirs('{}residual_b/'.format(plotDir),exist_ok=True)
		# plt.close('all')
		plt.figure(num=1,figsize=(6,6),clear=True)
		quiv_scale=100
		quiv_label_val = 5
		quiv_label = '{} pix'.format(quiv_label_val)
		q = plt.quiver(xts[include],yts[include],(xref[include]-xts[include]),(yref[include]-yts[include]),np.arctan2(yref[include]-yts[include], xref[include]-xts[include]),norm=Normalize(vmin=-math.pi,vmax=math.pi),  cmap='hsv', scale=quiv_scale, angles='xy',width=0.0005)
		# q = plt.quiver(xts[outliers],yts[outliers],(xref[outliers]-xts[outliers]),(yref[outliers]-yts[outliers]), color='red', scale=quiv_scale, angles='xy',width=0.0005)
		plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
		plt.xlim(-400,2448)
		plt.ylim(-400,2448)
		# plt.axis('equal')
		# plt.set_aspect('equal','box')
		plt.title('Distortion residuals_b {}'.format(ref_iteration))
		# plt.savefig('{}/residual_b/residual_b_quiver_{}_{}.pdf'.format(plotDir,solution_year,ref_iteration), bbox_inches='tight',dpi=200)
		plt.savefig('{}/residual_b/residual_b_quiver_{}_{}_{}.jpg'.format(plotDir,solution_year,ref_iteration,fp_iteration), bbox_inches='tight',dpi=600)


		for night in obs_nights:
		#plot individual frame residuals
			os.makedirs('{}resid_b_individual/{}/{}/{}'.format(plotDir,ref_iteration,night,fp_iteration),exist_ok=True)
			for f in range(len(osiris_filenames_dict[night])):
				# idx = np.where(frame[include] == f)[0]
				idx = np.where((frame[include] == f) & (night_c[include] == night))[0]
				plt.figure(num=1,figsize=(6,6),clear=True)
				quiv_scale=20
				quiv_label_val = 1
				quiv_label = '{} pix'.format(quiv_label_val)
				# print(f'idx = {idx}')
				angle_colour = np.arctan2(yref[include][idx]-yts[include][idx], xref[include][idx]-xts[include][idx])
				# norm = Normalize(vmin=0,vmax=2*math.pi)
				# norm.autoscale(angle_colour)
				# colourmap = 'hsv'
				q = plt.quiver(xts[include][idx],yts[include][idx],(xref[include][idx]-xts[include][idx]),(yref[include][idx]-yts[include][idx]),angle_colour, norm=Normalize(vmin=-math.pi,vmax=math.pi), cmap='hsv', scale=quiv_scale, angles='xy',width=0.002)
				# q = plt.quiver(xts[outliers],yts[outliers],(xref[outliers]-xts[outliers]),(yref[outliers]-yts[outliers]), color='red', scale=quiv_scale, angles='xy',width=0.0005)
				plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
				plt.xlim(-400,2448)
				plt.ylim(-400,2448)
				plt.xlabel('Pixels')
				plt.ylabel('Pixels')
				# plt.axis('equal')
				# plt.set_aspect('equal','box')
				plt.title('Distortion residuals_b {} {}'.format(ref_iteration, f))
				plt.savefig('{}resid_b_individual/{}/{}/{}/residual_b_quiver_{}_{}_{}.pdf'.format(plotDir,ref_iteration,night,fp_iteration,night,ref_iteration,f), bbox_inches='tight',dpi=200)


			#-----------------------------------

		# plt.close('all')

		xts_ra = np.array(xts_ra)		#distortion corrected observations transformed into RA Dec coordinates [include] (not centred)
		yts_dec = np.array(yts_dec)

		residuals_b_radec = np.hypot(raref[include]-xts_ra, decref[include]-yts_dec)
		median_4p_residual_radec = np.median(residuals_b_radec)

		residuals_b = np.hypot(xref-xts,yref-yts)     #distance from reference to transformed
		median_4p_residual = np.median(residuals_b[include])

		print(f'Ref_Iteration:{ref_iteration} fp_iteration:{fp_iteration} Median_residual:{median_4p_residual:.5f}')
		with open('{}fp_iteration_residuals_{}.txt'.format(resultDir,solution_year), 'a') as temp:
				temp.write(f'Ref_Iteration:{ref_iteration} fp_iteration:{fp_iteration} Median_residual_pix:{median_4p_residual:.5f} Median_residual_radec:{median_4p_residual_radec:.7f}\n')




	# return legendre_transformation ???
	#------------------------fp_iteration ends here---------------------------


	#need some parameters back here. current_distortion_solution and a bunch of parameters

	if centred:
		xref = xref + x_central 		# x_ref is already centred, so un-centre it.
		yref = yref + y_central			
		xc, yc = current_distortion_correction.evaluate(1024,1024)
		xc -= 1024
		yc -= 1024
	else:
		xc = 0
		yc = 0

	grid = np.arange(0,2048+1,64)
	
	xx, yy = np.meshgrid(grid, grid)
	x2,y2 = tform_leg.evaluate(xx,yy)
	x1 = xx.flatten()
	y1 = yy.flatten()
	x2 = x2.flatten()
	y2 = y2.flatten()

	dr = np.sqrt((x2-x1)**2 + (y2-y1)**2)
	print("Mean distortion before shift", np.mean(dr))

	a,b = tform_leg.evaluate(1024,1024)
	x_central = a-1024
	y_central = b-1024

	x2_c = x2 - x_central
	y2_c = y2 - y_central

	dr = np.sqrt((x2_c-x1)**2 + (y2_c-y1)**2)
	print("Mean distortion after shift", np.mean(dr))

	xts, yts = tform_leg.evaluate(x,y)

	distances2 = (xref-xts)**2 + (yref-yts)**2
	distances_ref = np.hypot(xref-x,yref-y)   #distance from OSIRIS original to reference.
	distances_tr = np.hypot(xts-x,yts-y)   #distance from OSIRIS original to transformed.
	unexplained_error = np.sqrt(unexplained_x_error**2 + unexplained_y_error**2)
	# unexplained_error = 0
	variances = 1/weights**2 + unexplained_error**2
	weights = 1/np.sqrt(variances)
	residuals_b_std = np.hypot((xts-xref),(yts-yref))*weights*0.01    #hoping that the units are wrong     10 mas per pixel.    0.01 arcsec per pixel

	#xe*0.01 to get it into arcseconds? 
	#* 0.01 to get it into pixels
	residuals_b_x = (xts-xref) / np.sqrt((xe*0.01)**2 + xeref**2 + unexplained_x_error**2) *0.01    #np.hypot(xe,xeref)
	residuals_b_y = (yts-yref) / np.sqrt((ye*0.01)**2 + yeref**2 + unexplained_y_error**2) *0.01    #np.hypot(ye,yeref)

	weighted_mean_residual_b = np.average(residuals_b[include], weights = weights[include])
	weighted_mean_residual_b_squared = np.sqrt(np.average(residuals_b[include]**2, weights = weights[include]))

	print('min,median,max residuals:',np.min(residuals_b[include]),np.median(residuals_b[include]),np.max(residuals_b[include]))
	print('min,median,max uncertainties:',np.min(1/weights[include]),np.median(1/weights[include]),np.max(1/weights[include]))
	print('weighted mean residual_b: {}'.format(weighted_mean_residual_b))
	print('weighted mean residual_b squared: {}'.format(weighted_mean_residual_b_squared))

	bad_a = np.nonzero(residuals_b > 22)
	bad_o = np.nonzero(outliers)
	# print(bad_a)
	# print(bad_o)																		#manually finding outliers by their residuals
	print('Residual_b > 22, not caught as outliers:', np.setdiff1d(bad_a, bad_o, assume_unique=True))
	# quit()

	#----------------------------------

	# median_residual_b = np.median(residuals_b[include])
	# median_residuals_b.append(median_residual_b)
	median_residuals_b.append(median_4p_residual)
	mean_residuals_b.append(weighted_mean_residual_b)
	mean_residuals_b_squared.append(weighted_mean_residual_b_squared)

	min_residuals_b.append(np.min(residuals_b[include]))
	max_residuals_b.append(np.max(residuals_b[include]))
	num_residuals_b.append(len(residuals_b[include]))

	median_residuals_b_radec.append(median_4p_residual_radec)

	#This section was at the end, after A. Moved here.
	with open(resultDir + 'iteration_residuals_{}.txt'.format(solution_year), 'w') as temp:
		for r, resid in enumerate(median_residuals_b):
			# temp.write(str(resid) + '\n')
			temp.write(f'Min:{min_residuals_b[r]:.5f}  Median:{median_residuals_b[r]:.5f}  Max:{max_residuals_b[r]:.5f}  Num:{num_residuals_b[r]:.5f}  Weighted mean:{mean_residuals_b[r]:.5f}  Weighted mean squared:{mean_residuals_b_squared[r]:.5f} Median_mas:{median_residuals_b_radec[r]:.7f} \n') #| Mean_a:{mean_residuals_a[r]:.5f}

	print('Median residual = {}'.format(median_4p_residual))
	print('Median residual radec = {}'.format(median_4p_residual_radec))		
	print('Weighted mean residual = {}'.format(weighted_mean_residual_b))
	print('Weighted mean residual squared = {}'.format(weighted_mean_residual_b_squared))
	#removed break condition. Will complete all ref_iterations
	# if ref_iteration > 5:
	# 	# improvement = median_residuals_b[ref_iteration-1]-median_residuals_b[ref_iteration]
	# 	improvement = mean_residuals_b[ref_iteration-1]-mean_residuals_b[ref_iteration]
	# 	print('Improvement = {}'.format(improvement))
	# 	if ref_iteration > 5:
	# 		if 0 < improvement < 0.005:					#0.05 mas =  0.005 pixels    50 mas = 5 pixel. 0.05 mas = 0.005
	# 			print('Improvement < 0.05, breaking at ref_iteration {}'.format(ref_iteration))
	# 			break
	print('All ref_iteration median residuals_b = {}'.format(median_residuals_b))
	print('All ref_iteration weighted mean residuals_b = {}'.format(mean_residuals_b))
	print('All ref_iteration weighted mean residuals_b squared= {}'.format(mean_residuals_b_squared))
	print('All ref_iteration mean residuals_a = {}'.format(mean_residuals_a))
	# plt.close('all')
	# quit()

def load_osiris_file(starfindDir,filename):
	#Loads a starlist file, converts it to a StarList object.
	with open (starfindDir + filename) as starfindFile:
			starfindTable = starfindFile.read().splitlines()
	osirisStarList = starlists.StarList.from_lis_file(starfindTable, error=False)  #the starlists may be in _stars.txt files, not .lis files.
	osirisStarList = starlists.StarList.from_table(starfindTable, error=False) 

	#from_lis_file requires these columns in order: name, m, t, x ,y , snr, corr, N_frames, flux
	# these are in the _stars.txt files, except for snr and corr.
	# modify my other script to generate appropriate files.
	#from_table requires name, x, y, and m.

	return osirisStarList


def reference_section_plots(refTable_b):
	plt.figure(num=1,figsize=(6,6),clear=True)
	ptsize = 5
	os.makedirs('{}combined_ref'.format(plotDir),exist_ok=True)
	# xbound = np.logical_and(-4.99 < refTable_b['x0'],refTable_b['x0'] < -4.77)
	# ybound = np.logical_and(8.17 < refTable_b['y0'],refTable_b['y0'] < 8.4)
	# s = np.where(np.logical_and(xbound,ybound))[0]
	star_id = '86483'
	s = np.where(refTable_b['name'] == star_id)[0]
	plt.scatter(refTable_b['x0'][s],refTable_b['y0'][s],ptsize,'k',label='Mean')
	plt.scatter(refTable_a['x'][s],refTable_a['y'][s],ptsize,label='Individual frames')
	# plt.scatter(refTable_a['x'][s][0][0],refTable_a['y'][s][0][0],ptsize,label='Hubble?')
	
	s = np.where(refTable_current[night]['name'].astype('str') == star_id)[0]
	plt.scatter(refTable_current[night]['x0'][s],refTable_current[night]['y0'][s],ptsize,label='Ref')
	plt.xlim([-4.932,-4.88])
	plt.ylim([8.18,8.242])
	plt.title('Averaging star position it:{}'.format(ref_iteration))
	plt.legend(loc='upper left')
	plt.xlabel('Relative RA (arcsec)')
	plt.ylabel('Relative Dec (arcsec)')
	plt.savefig('{}combined_ref/new_reference_stars_{}_{}.jpg'.format(plotDir,solution_year,ref_iteration), bbox_inches='tight',dpi=200)
	# plt.show()

	print('{},{} matched stars in refTable_a'.format(len(idx1),len(idx2)))
	#plot individual frame residuals
	os.makedirs('{}residual_A/'.format(plotDir,ref_iteration),exist_ok=True)
	# for f in range(len(osiris_filenames)):
	plt.figure(num=1,figsize=(6,6),clear=True)
	# plt.clf()
	quiv_scale=0.2
	quiv_label_val = 0.01
	quiv_label = '{} arcsec'.format(quiv_label_val)
	# idx = np.where(frame[include] == f)[0]
	# print(f'idx = {idx}')
	# q = plt.quiver(xts[include][idx],yts[include][idx],(xref[include][idx]-xts[include][idx]),(yref[include][idx]-yts[include][idx]),np.arctan2(yref[include][idx]-yts[include][idx], xref[include][idx]-xts[include][idx]),norm=Normalize(vmin=-math.pi,vmax=math.pi),  cmap='hsv', scale=quiv_scale, angles='xy',width=0.005)
	# q = plt.quiver(xts[outliers],yts[outliers],(xref[outliers]-xts[outliers]),(yref[outliers]-yts[outliers]), color='red', scale=quiv_scale, angles='xy',width=0.0005)
	q = plt.quiver(refTable_b['x0'][idx1],refTable_b['y0'][idx1],(refTable_current[night]['x0'][idx2]-refTable_b['x0'][idx1]),(refTable_current[night]['y0'][idx2]-refTable_b['y0'][idx1]),np.arctan2(refTable_current[night]['y0'][idx2]-refTable_b['y0'][idx1], refTable_current[night]['x0'][idx2]-refTable_b['x0'][idx1]), norm=Normalize(vmin=-math.pi,vmax=math.pi), cmap='hsv', scale=quiv_scale, angles='xy',width=0.005)
	plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
	plt.xlim(-20,20)
	plt.ylim(-20,20)
	plt.xlabel('relative RA (arcsec)')
	plt.ylabel('relative Dec (arcsec)')
	# plt.axis('equal')
	# plt.set_aspect('equal','box')
	plt.title('Distortion residuals_a it:{}'.format(ref_iteration))
	plt.savefig('{}residual_A/residual_quiver_{}_{}.jpg'.format(plotDir,solution_year,ref_iteration), bbox_inches='tight',dpi=200)



def get_PA(filename):
	hdr = fits.getheader(filename)
	# return math.radians(hdr['PA_IMAG'])
	return hdr['PA_IMAG']

def calculate_PCU_guess(filename):
	if config['instrument'] == 'OSIRIS':
		hdr = fits.getheader(filename)
		# pcu_x = float(hdr['PCSFX'])   #these keywords may change
		# pcu_y = float(hdr['PCSFY'])
		# pcu_r = hdr['PCUR']
		# pcu_r = 65.703
		raw_filename = filename[-25:-10] + filename[-5:]
		mask = log['Filename'] == raw_filename
		pcu_params = log[mask]

		pcu_x = pcu_params['x']
		pcu_y = pcu_params['y']
		pcu_z = pcu_params['z']
		pcu_r = pcu_params['r']
		# print(pcu_x,pcu_y,pcu_r)

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

		S = 138.5  #pixels per mm
		theta_offset = 65.703   #this is some offset in degrees between the reported PCU rotation and the orientation we want.
		theta = np.radians(pcu_r-theta_offset)  #check the sign of the rotation
		a1 = S * math.cos(theta)
		a2 = S * math.sin(theta)

		pinhole_scale = 138.5 #pixles per mm on the pinhole mask.
		pinhole_x_offset = 91.05  # mm
		pinhole_y_offset = 183.95  # mm
		pinhole_x_center = 1667 +9 #pixels.  The location of the center of the pinohole mask in pixels when the PCU is at the offset position above.
		# pinhole_y_center = 613 #pixels  On raw image
		pinhole_y_center = 1422 #on flipped image

		a0 = S*(pcu_x - pinhole_x_offset) + pinhole_x_center  #in pixels
		b0 = S*(pinhole_y_offset - pcu_y) + pinhole_y_center #flip applied to image, so using negative pix coordinates

		four_p = Empty()
		four_p.__class__ = transforms.four_paramNW   #I am trying to create an instance of the class without calling __init__, which requires two lists of stars to be passed in.
		four_p.px = [a0,a1,a2]
		four_p.py = [b0,-a2,a1]
		four_p.order = None
		four_p.mag_offset = 0
	else:
		print('Instruments other than OSIRIS not supported')
		sys.exit(0)
	return four_p

class Empty(object):
	pass


def mag_cut(starlist,min,max):
	idx = np.where((starlist['m'] >= min) & (starlist['m'] <= max))
	return starlist[idx]


def stack_frames():
	pass



def load_reference_list():
	#remember, there can be different lists in one run
	pass


def load_observed_lists():
	pass



main()









