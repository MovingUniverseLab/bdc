#distortion.py

#this script  reads in the stack of star positions and the gaia data, and does the thing.
# import sys
# print('\n'.join(sys.path))

import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table
# from astroquery.gaia import Gaia
import pickle
import os
import fnmatch
from astropy.io import ascii, fits
from astropy.stats import sigma_clipping
from astroquery.gaia import Gaia
from astroquery.mast import Observations, Catalogs
# from flystar import align, match, transforms
# import flystar
from flystar import analysis, align, transforms, starlists, plots
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg') #stops memory leak?
from matplotlib.colors import LogNorm, Normalize
import numpy as np
import math
# import AstroAtmosphere as AAt 
from kai.reduce import dar
from kai import instruments
import datetime
# from jlu.observe import weather
from collections import Counter
import requests
# import memory_profiler as mp
from memory_profiler import profile 	#do not import when running timed mode.

# from header of ci200804_a014002_flip.fits
# TARGRA  =           322.492875 / target right ascension (21:29:58.29 h)         
# TARGDEC =          12.16679167 / target declination (+12:10:00.4 deg)
# RA      =         322.48999069 / telescope right ascension (21:29:57.60 h)      
# DEC     =          12.17453385 / telescope declination (+12:10:28.3 deg)        

# from header of ci200804_a015002_flip.fits
# RA      =         322.48989305 / telescope right ascension (21:29:57.57 h)      
# DEC     =           12.1758901 / telescope declination (+12:10:33.2 deg)        


# PSCALE  =                 0.01 / Pixel scale for science camera    (arcsec, = 10 mas)             
# NAXIS1  =                 2048 / Axis length                                    
# So FoV = 20.48 arcseconds
#0.0015 deg = 5.4 arcsec. So total fov = 15 arcsec from centre. = 30x30 = 0.00833 deg

# s1 ra 322.48982502490475        dec 12.17430185928072
osiris_scale = 0.009942  #arcseconds per pixel


#---------------------------------------------
recalculate_plate_scale = False   #re-calculate the 4-parameter fit with the distortion corrected data to calculate the plate scale.
create_stacks = False

# run_label = 'used_in_trans'   
# run_label = 'order' #moved to inside main()

create_combined_reflist = True		#if True, will force recalculation of the combined reflist (section A). If false, will attempt to load from file first
fit_to_combined_reflist = True 		#if True, will force recalculation of the fit to the reflist (section B). If false, will attempt to load from file first
previous_results_location = './run_3/results_stacked/'  #location to search for previous combined reflist and fit.
# transform_files_location = './transform_files/'	#location to search for previous transformation guesses.
# transform_files_location = './initial_transformations/Hubble_FF/'	#testing an older version of the transform files to see if it helps.
transform_files_location = './run_3/transformations/Hubble_FF/'	

# reference_instrument = 'GAIA'
reference_instrument = 'Hubble'
#---------------------------------------------

centred = True
debugging = False
use_flystar_velocity = False
manually_average_star_positions = True
trim_not_used_in_trans = False

# @profile
def main(solution_year,fitmode='FF',order=5):
	print('Starting solution for', solution_year)
	run_label = 'full_1_order_' + str(order)

	# global nightDir, cleanDir, use_prev_gaia, targetID, ra_field, dec_field,radius, minmag, single_fit, rad_tolerance 
	# global mag_tolerance, mag_limits, bad_files, dont_trim, sl, show_plots, starfindDir, plotDir, plotDir_n, resultDir, stackDir, tformDir, target

	global nightDir, cleanDir, use_prev_gaia, targetID, ra_field, dec_field,radius, minmag, single_fit, rad_tolerance 
	global mag_tolerance, mag_limits, bad_files, dont_trim, sl, show_plots, starfindDir, plotDir, plotDir_n, resultDir, stackDir, tformDir

	if solution_year == '2020':
		obs_nights = ['n1','n2','n3','n4']
	elif solution_year == '2021':
		obs_nights = ['n5']
	else:
		print('Year must be 2020 or 2021')
		quit()


	target_dict = {
	'n1':'m15',
	'n2':'m92',
	'n3':'m92',
	'n4':'m92',
	'n5':'m15',
	}
	nightDir_dict = {
	'n1':'/u/mfreeman/work/d/n1/',
	'n2':'/u/mfreeman/work/d/n2/',
	'n3':'/u/mfreeman/work/d/n3/',
	'n4':'/u/mfreeman/work/d/n4/',
	'n5':'/u/mfreeman/work/d/n5/',
	}
	cleanDir_dict = {
	'n1':nightDir_dict['n1'] + 'clean/m15_kn3_tdOpen/',
	'n2':nightDir_dict['n2'] + 'clean/m92_kp_tdOpen/',
	'n3':nightDir_dict['n3'] + 'clean/m92_kp_tdOpen/',
	'n4':nightDir_dict['n4'] + 'clean/m92_kp_tdOpen/',
	'n5':nightDir_dict['n5'] + 'clean/m15_kn3_tdhBand/',
	}

	starfindDir_dict = {
	'n1':cleanDir_dict['n1'] + 'starfinder/',
	'n2':cleanDir_dict['n2'] + 'starfinder/',
	'n3':cleanDir_dict['n3'] + 'starfinder/',
	'n4':cleanDir_dict['n4'] + 'starfinder/',
	'n5':cleanDir_dict['n5'] + 'starfinder/',
	}
	stackDir_dict = {
	'n1':cleanDir_dict['n1'] + 'stacks/',
	'n2':cleanDir_dict['n2'] + 'stacks/',
	'n3':cleanDir_dict['n3'] + 'stacks/',
	'n4':cleanDir_dict['n4'] + 'stacks/',
	'n5':cleanDir_dict['n5'] + 'stacks/',
	}

	osiris_filenames_dict = {
	'n1':get_osiris_files(stackDir_dict['n1']),
	'n2':get_osiris_files(stackDir_dict['n2']),
	'n3':get_osiris_files(stackDir_dict['n3']),
	'n4':get_osiris_files(stackDir_dict['n4']),
	'n5':get_osiris_files(stackDir_dict['n5']),
	}

	bad_files_dict = {
	'n1':['ci200804_a022007_flip_0.8_stf.lis','ci200804_a026012_flip_0.8_stf.lis','ci200804_a027003_flip_0.8_stf.lis',],
	'n2':[],
	'n3':[],
	'n4':['ci200814_a032002_flip_0.8_stf.lis','ci200814_a032007_flip_0.8_stf.lis','ci200814_a032008_flip_0.8_stf.lis'],  #the first one missed a bright star, so messes up the matching for stacking.
	'n5':['ci211024_a011009_flip_0.8_stf.lis','ci211024_a012010_flip_0.8_stf.lis','ci211024_a020012_flip_0.8_stf.lis'],
	}

	mag_limits_dict = {
	'n1':[6,16],
	'n2':None,
	'n3':None,
	'n4':[6,16],
	'n5':[6,16],
	}

	mag_tolerance_dict = {
	'n1':[2, 2, 2],
	'n2':[2, 2, 2],
	'n3':[2, 2, 2],
	'n4':[2, 2, 2],
	'n5':[2, 2, 2], 
	}

	plotDir = './' + run_label + '/plots_stacked/' + reference_instrument + '_' + fitmode+ '/'
	resultDir = './' + run_label + '/results_stacked/' + reference_instrument + '_' + fitmode + '/'
	tformDir = './' + run_label + '/transformations/' + reference_instrument + '_' + fitmode + '/'

	os.makedirs(plotDir, exist_ok=True)
	os.makedirs(resultDir, exist_ok=True)
	os.makedirs(tformDir, exist_ok=True)


	bad_gaia_stars = []
	# bad_gaia_stars = [1745948328025672960, 1745948328028860800]#, 1745948396745203456, 1745948396746376448, 1745948396748336384, 1745948396750867072]
						#1745948323734090368, is the target, so can't flag as bad
	bad_gaia_stars = [1360405469097523456, 1360405469102049792, 1360405503457185792, 1360405503457537408, 1360405503461809536, 1360405503461940736, 
						1360405503461949056, 1360405503461967744, 1745948328025673472, 1745948328025674240, 1745948328025687296, 1745948328028724736, 
						1745948396745149312, 1745948396745157120, 1745948328025321472]


	# --------------------------------------------------------------------




	if create_stacks:  
		print('Creating new stacked images')
		for night in obs_nights:
		# for night in ['n3']: 
			print('Night', night)
			previous_coo = np.array([0,0])
			stack_list = {}
			stacknum = 0

			nightDir = nightDir_dict[night]	
			cleanDir = cleanDir_dict[night]
			starfindDir = cleanDir + 'starfinder/'
			stackDir = cleanDir + 'stacks/'
			osiris_filenames = get_osiris_files(starfindDir)
			# osiris_filenames = osiris_filenames[sl]
			# osiris_filenames = ['ci211024_a009002_flip_0.8_stf.lis','ci211024_a009003_flip_0.8_stf.lis','ci211024_a009004_flip_0.8_stf.lis']
			print(len(osiris_filenames), 'OSIRIS images to stack')

			for i, filename in enumerate(osiris_filenames):
				if filename in bad_files_dict[night]:
					print('{} {} flagged as bad, skipping'.format(i,filename))
					continue
					# errorlist.append(filename[:-4])
				else:
					print(filename)
					# starlist = load_osiris_file(filename)
					# plt.close('all')

					with open (cleanDir + filename[:-12] + '.coo') as coofile:
						cooxy = np.genfromtxt(coofile)

					coo_shift = previous_coo - cooxy
					if np.hypot(coo_shift[0],coo_shift[1]) > 5:
						#start a new stack
						stackname = str(stacknum)
						stacknum += 1
						stack_list[stackname] = [filename]
					else:
						#append last stack
						stack_list[stackname].append(filename)
						#add to current stack
					previous_coo = cooxy



			trans_guess_location = './stack_guesses/'#.format(night)
			os.makedirs(trans_guess_location, exist_ok=True)


			for i, stack in enumerate(stack_list):
				nframes = len(stack_list[stack])
				print(i, 'has', nframes,'frames:', stack_list[stack])
				# if i < 36: 
				# 	continue
				if nframes < 2:
					print('Only', nframes, 'images, skipping')
				else:
					
					list_of_starlists_s = [load_osiris_file(starfindDir ,filename) for filename in stack_list[stack]]

					print('Number of stars in each list:', [len(j) for j in list_of_starlists_s])
					# for i, starlist in enumerate(list_of_starlists_s):
					# 	print(stack_list[stack][i], len(starlist))

					trans_guess_files = []
					trans_guesses = []
					for j in range(nframes):
						trans_guess_files.append('{}_s{}_{}'.format(night,i,j))
					
					# if all([os.path.exists(trans_guess_location+trans_guess_files[j]) for j in range(nframes)]):
					# 	for j in range(nframes):
					# 		print('Loading transform guess', trans_guess_files[j])						
					# 		with open(trans_guess_location+trans_guess_files[j], 'rb') as temp:
					# 			trans_guesses.append(pickle.load(temp)) 					#initial guess for the transformation
					# else:
					# 	# print('Not loading transform guess')
					# 	# trans_guesses= None
					
					for j in range(nframes):
						trans_guesses.append(transforms.PolyTransform(1,[0,1,0],[0,0,1]))


					#dr_tol 0.4
					msc = align.MosaicSelfRef(
						list_of_starlists_s, iters=3,
						dr_tol=[2.0,2.0,2.0], dm_tol=mag_tolerance_dict[night],
						outlier_tol=[None, None, None],
						trans_class=transforms.PolyTransform,
						trans_input=trans_guesses,
						# trans_class=transforms.four_paramNW,
						trans_args=[{'order': 0}, {'order': 0}, {'order': 0}],
						use_vel=False,
						mag_trans=True,
						mag_lim=mag_limits_dict[night], #[6,16],
						weights=None,
						calc_trans_inverse=True,    
						init_guess_mode='miracle', verbose=0)

					msc.fit()
					tab1 = msc.ref_table
					# tform = msc.trans_list
					# tform_inv = msc.trans_list_inverse 
					for j in range(nframes):
						with open(trans_guess_location+trans_guess_files[j], 'wb') as temp:
							pickle.dump(msc.trans_list[j], temp)
			

					idx = np.where(np.all(tab1['N_frames'] == 1,axis=1,))[0]  #only select stars that were seen in all frames.
					print(len(idx), 'stars found in all frames')
					# idx = np.where(np.sum(tab1['N_frames'][:],axis=1)> 2)[0]

					idx = np.where(np.count_nonzero(tab1['N_frames']==1,axis=1) > nframes/2)[0]
					print('{} stars found in >{} frames'.format(len(idx),nframes/2))
					# for o , thing in enumerate(tab1['N_frames']):
					# 	print(thing)
					# print(tab1['N_frames'][:,:],tab1['N_frames'][0][1],tab1['N_frames'][0][2])

					stacked_starlist = tab1['name',][idx]
					stacked_starlist['x'] = tab1['x0'][idx]
					stacked_starlist['xe'] = tab1['x0e'][idx]
					stacked_starlist['y'] = tab1['y0'][idx]
					stacked_starlist['ye'] = tab1['y0e'][idx]
					stacked_starlist['m'] = tab1['m0'][idx]
					stacked_starlist['me'] = tab1['m0e'][idx]

					# stacked_starlist['t'] = tab1['t'][:,-1][idx]
					# stacked_starlist['snr'] = tab1['snr'][:,-1][idx]
					# stacked_starlist['corr'] = tab1['corr'][:,-1][idx]
					# stacked_starlist['flux'] = tab1['flux'][:,-1][idx]
					stacked_starlist['t'] = np.nanmean(tab1['t'],axis=1)[idx]
					stacked_starlist['snr'] = np.nanmean(tab1['snr'],axis=1)[idx]
					stacked_starlist['corr'] = np.nanmean(tab1['corr'],axis=1)[idx]
					stacked_starlist['flux'] = np.nanmean(tab1['flux'],axis=1)[idx]

					# stacked_starlist['N_frames'] = np.sum(tab1['N_frames'][:],axis=1)[idx]
					stacked_starlist['N_frames'] = np.count_nonzero(tab1['N_frames']==1,axis=1)[idx]

					starlists.StarList.to_lis_file(stacked_starlist,stackDir + stack_list[stack][0] )


		#night loop ends here
		print('Stacked creation complete')
		quit()

 	#-----------------------------------------Finished creating stacks------------------------------------------

	if reference_instrument == 'Hubble':

		hubble_files = {'m15':'/g/lu/data/m15/hst_ref/NGC7078cb.pm','m92':'/g/lu/data/m92/hst_ref/NGC6341cp.pm'}

		if solution_year == '2021':
			hubbleData = {'m15':fetch_hubble(hubble_files['m15'])}
		elif solution_year == '2020':
			hubbleData = {'m15':fetch_hubble(hubble_files['m15']), 'm92':fetch_hubble(hubble_files['m92'])}


	elif reference_instrument == 'GAIA':
		#out of date
		gaiaData = fetch_gaia(ra_field, dec_field, radius, targetID, night, use_prev_gaia)
		
		mask = np.where(np.isin(gaiaData['source_id'], bad_gaia_stars))[0]
		gaiaData.remove_rows(mask)
		print("Removed", np.count_nonzero(mask), "Gaia stars that I flagged as bad")

		max_astrometric_excess_noise = 2
		mask = gaiaData['astrometric_excess_noise'] > max_astrometric_excess_noise
		gaiaData.remove_rows(mask)
		print("Removed", np.count_nonzero(mask), "Gaia stars with >", max_astrometric_excess_noise, "mas of excess astrometric noise")
		print(len(gaiaData), 'Gaia stars remaining')


		
		starlist0 = load_osiris_file(stackDir ,osiris_filenames[0])
		gaiaData = project_pos(gaiaData,starlist0,'gaia')

		#instead of projecting the positions myself, set the Flystar flag use_vel=True?
		#No, I apply DAR before running flystar, so I need the positions corrected first.
	 
		cooi = np.where(gaiaData['source_id'] == targetID)[0]
		refTable_H = analysis.prepare_gaia_for_flystar(gaiaData, gaiaData['ra'][cooi], gaiaData['dec'][cooi])

	else:
		print(reference_instrument, 'should be Hubble or GAIA')
		quit()


    

	if debugging:
		ref_iterations = [0,1,2,3]
	else:
		ref_iterations = [0,1,2,3,4,5,6,7]

		# ref_iterations = [0]
		# ref_iterations = [0,1,2,3]
	# ref_iteration_length = 4
	correction_list = []	
	
	# for i in range(ref_iteration_length):
	# 	correction_list.append(load_legendre_iteration(solution_year,i))

	# xc, yc = correction_list[0].evaluate(1024,1024)
	# xc -= 1024
	# yc -= 1024


	fp_iterations = [0,1,2,3,4,5]
	# fp_iterations = [0,1,2]

	osiris = instruments.OSIRIS()
	

	# open('{}fp_iteration_residuals_{}.txt'.format(resultDir,solution_year), 'w').close()	 #clear the file with the 4p loop residuals.
	with open('{}fp_iteration_residuals_{}.txt'.format(resultDir,solution_year), 'w') as temp:
		pass

	median_residuals_b = []
	mean_residuals_b = []
	mean_residuals_b_squared = []
	min_residuals_b = []
	max_residuals_b = []
	num_residuals_b = []

	median_residuals_b_radec = []

	mean_residuals_a = []
	
	refTable_current = {}	#holds most recent refTable for each night
	refTable_current_filename = '{}refTable_current_{}'.format(resultDir,solution_year) 

	#-----------------------------Star reference iteration here-------------------------------------------
	for ref_iteration in ref_iterations:
		print('\n \n Ref Iteration {} \n'.format(ref_iteration))
		# plt.close('all')









		# --------------------------------------- Section B ---------------------------------------------------------------


		#fp_iteration starts here?
		#get the distortion solution, correct the data, pass that back into the fitter.
		#Do I want to save all fp_iterations, or just save the final fit? I think I want all iterations to be saved. So I can quickly re-run and plot them.


		current_distortion_correction = None
		tab1_initial = {'n1':[],'n2':[],'n3':[],'n4':[],'n5':[]} #is regenerated each ref_iteration
		#if I am loading dist files then this can break if a set of fp_iterations was incomplete - it won't run the first iteration to generate it. Need to save as a file?

		for fp_iteration in fp_iterations:

			# plt.close('all')


			print(f'\n Ref Iteration={ref_iteration} 4p Iteration={fp_iteration} \n')

			dist_file = '{}dist_measures_{}_{}_{}.txt'.format(resultDir,solution_year,ref_iteration,fp_iteration)
			dist_file_previous = '{}dist_measures_{}_{}_{}.txt'.format(resultDir,ref_iteration,fp_iteration,solution_year) #old name convention
			# dist_file_previous = '{}dist_measures_{}_{}_{}.txt'.format(previous_results_location + reference_instrument + '_' + fitmode + '/',ref_iteration,fp_iteration,solution_year) 

			if os.path.exists(dist_file):
				dist_file_toUse = dist_file
			elif os.path.exists(dist_file_previous):
				dist_file_toUse = dist_file_previous
			else:
				dist_file_toUse = None

			# if (not fit_to_combined_reflist) and os.path.exists(dist_file):
			# 	print('b: Loading fit {}'.format(dist_file))
			# 	distortion_data = ascii.read(dist_file,format='fixed_width')


			if (not fit_to_combined_reflist) and (dist_file_toUse is not None) and (os.path.exists(refTable_current_filename)):
				print('b: Loading fit {}'.format(dist_file_toUse))
				distortion_data = ascii.read(dist_file_toUse,format='fixed_width')


			else:
				print('b: Generating new fit')
				errorlist = []
				successlist = []
				x_O = []
				y_O = []
				x_G = []
				y_G = []
				weights = []
				idnums = []
				starlist_num = []
				xe_O = []
				ye_O = []
				xe_G = []
				ye_G = []
				used_in_trans_1 = []
				Ra_O = []
				Dec_O = []
				Ra_G = []
				Dec_G = []
				night_col = []


				used_files = []
				scales = []
				rotations = []
				PAs = []

				#-------------------------Night loop starts here-------------------------------
				for night in obs_nights:
					print('Starting night', night)
					if night == 'n1':
						hubble_file = '/g/lu/data/m15/hst_ref/NGC7078cb.pm' #'je0o61lzq_flt.xymrduvqpk'
						targetID = 1745948323734090368   #gaia ID of target star. Ra and Dec used in prepare_gaia_for_flystar()
						# ra_field = '21:29:57.60'  #approximate centre of FoV, selects gaia stars in radius  
						# dec_field = '12:10:28.3'
						ra_field = 322.48999069
						dec_field = 12.17453385
						radius = 20   #arcseconds
						minmag = 15.4  #dimmest mag for cut
						single_fit = True  #run Mosaic2Ref for each image individually
						rad_tolerance = [0.4, 0.4, 0.2]
						dont_trim = ['ci200804_a014004_flip_0.8_stf.lis', 'ci200804_a026009_flip_0.8_stf.lis']
						# 026002 is bad either way. Didn't find the brightest star. Could be that flag I selected. Try without?
						sl = slice(0,None)
						# sl = slice(34,38)
						show_plots = False
					elif night == 'n2':
						hubble_file = '/g/lu/data/m92/hst_ref/NGC6341cp.pm' #'idk901xpq_flt.xymrduvqpk'
						targetID = 1360405503461790848  #gaia ID of target star. Ra and Dec used in prepare_gaia_for_flystar()
						ra_field = 259.285096306648     #approximate centre of FoV, selects gaia stars in radius
						dec_field = 43.13751895071527
						radius = 20   #arcseconds
						minmag = 15.4  #dimmest mag for cut
						rad_tolerance = [0.4, 0.4, 0.2]
						single_fit = True  #run Mosaic2Ref for each image individually
						dont_trim = []
						sl = slice(0,None)
						# sl = slice(151,None)
						show_plots = False    
					elif night == 'n3':
						hubble_file = '/g/lu/data/m92/hst_ref/NGC6341cp.pm' #'idk901xpq_flt.xymrduvqpk'
						targetID = 1360405503461790848  #gaia ID of target star. Ra and Dec used in prepare_gaia_for_flystar()
						ra_field = 259.285096306648     #approximate centre of FoV, selects gaia stars in radius
						dec_field = 43.13751895071527
						radius = 20   #arcseconds
						minmag = 15.4  #dimmest mag for cut
						rad_tolerance = [0.4, 0.4, 0.2]
						single_fit = True  #run Mosaic2Ref for each image individually
						dont_trim = []
						sl = slice(0,None)
						# sl = slice(151,None)
						show_plots = False    
					elif night == 'n4':
						hubble_file = '/g/lu/data/m92/hst_ref/NGC6341cp.pm' #'idk901xpq_flt.xymrduvqpk'
						targetID = 1360405503461790848  #gaia ID of target star. Ra and Dec used in prepare_gaia_for_flystar()
						ra_field = 259.285096306648     #approximate centre of FoV, selects gaia stars in radius
						dec_field = 43.13751895071527
						radius = 20   #arcseconds
						minmag = 15.4  #dimmest mag for cut
						rad_tolerance = [0.4, 0.4, 0.2]
						single_fit = True  #run Mosaic2Ref for each image individually
						dont_trim = []
						sl = slice(0,None)
						# sl = slice(151,None)
						show_plots = False  

					elif night == 'n5':
						hubble_file = '/g/lu/data/m15/hst_ref/NGC7078cb.pm' #'icbe05m9q_flt.xymrduvqpk'     #there are other files too?
						targetID = 1745948328028761984  #gaia ID of target star. Ra and Dec used in prepare_gaia_for_flystar()
						ra_field = 322.4912419147502  #approximate centre of FoV, selects gaia stars in radius
						dec_field = 12.164721331771977	
						radius = 20   #arcseconds
						minmag = 15.4  #dimmest starlist mag for cut
						rad_tolerance = [0.4, 0.4, 0.2]
						mag_tolerance = [2, 2, 2]
						mag_limits = [6,16]
						single_fit = True  #run Mosaic2Ref for each image individually
						dont_trim = []
						sl = slice(0,None)
						# sl = slice(79,None)
						show_plots = False 
					else:
						print('No night selected')
						quit()


					target = target_dict[night]
					nightDir = nightDir_dict[night]
					cleanDir = cleanDir_dict[night]
					starfindDir = starfindDir_dict[night]
					stackDir = stackDir_dict[night]
					osiris_filenames = osiris_filenames_dict[night]
					bad_files = bad_files_dict[night]
					mag_limits = mag_limits_dict[night]

					starfindDir = cleanDir + 'starfinder/'
					stackDir = cleanDir + 'stacks/'
					plotDir_n = plotDir + night + '/'

					make_dir(plotDir,ref_iteration)

					osiris_filenames = get_osiris_files(stackDir)

					sl2 = slice(0,None) 
					osiris_filenames = osiris_filenames[sl2]
					print(osiris_filenames)
					print(len(osiris_filenames), 'OSIRIS images')


					
					if night not in refTable_current.keys():

						if not use_flystar_velocity:
							starlist0 = load_osiris_file(stackDir ,osiris_filenames[0])
							hubbleData_p = project_pos(hubbleData[target],starlist0,'hubble_'+target)
							refTable_H = prepare_hubble_for_flystar(hubbleData_p,ra_field,dec_field,target)

						else:
							refTable_H = prepare_hubble_for_flystar(hubbleData,ra_field,dec_field,target)
						#generate osiris_filenames here

						# refTable_current = refTable_H.filled()		#unmask, required for quiver plots
					
						refTable_current[night] = refTable_H.filled()

						with open(refTable_current_filename, 'wb') as temp:
							pickle.dump(refTable_current, temp)


					for i, filename in enumerate(osiris_filenames):
						if filename in bad_files:
							print('{} {} flagged as bad, skipping'.format(i,filename))
							# errorlist.append(filename[:-4])
						else:

							starlist = load_osiris_file(stackDir ,filename)
							# plt.close('all')

							# ido = np.where(starlist['m'] < 15.5)
							# starlist = starlist[ido]
							fitsfile = cleanDir + filename[:-12] + '.fits'

							PA = get_PA(fitsfile)
							# print('PA', PA)
							# starlist = brightest_n(starlist,170)
							starlist = mag_cut(starlist,0,minmag)
							if not filename in dont_trim:
								starlist = edge_cut(starlist,5)

							if len(starlist) == 0:
								print(i,filename, '0 stars remaining after edge cut, skipping image')
								errorlist.append(filename[:-4])
								continue

							#----------------------------------------
							# plt.figure()
							# plt.hist(starlist['m'])
							# # plt.scatter(starlist['m'],starlist['vxe'],ptsize,alpha=0.2,label='Ref')
							# plt.show()
							# quit()
							#--------------------------------------

							if recalculate_plate_scale:
								xt, yt = correction.evaluate(starlist['x'],starlist['y'])
								# xt = xt.round(5)
								# yt = yt.round(5)
		

								# print(starlist['y'])
								starlist['x'] = xt 
								starlist['y'] = yt 



							# if reference_instrument == 'Hubble':
							# 	refTable_t = refTable
							# elif reference_instrument == 'GAIA':
							
							#I have use_vel = True when generating the combined reference frame, so the epoch is set to the osiris epoch, so no need to project positions.

							refTable_t = trim_gaia(refTable_current[night],filename,PA)    
							refTable_tf = refTable_t.filled()		#unmask, required for quiver plots

							plt.figure(num=4,figsize=(6,6),clear=True)

							if recalculate_plate_scale:
								refTable_d = dar.applyDAR(fitsfile, refTable_tf, plot=False, instrument=osiris, plotdir=plotDir_n + 'dar_a/'+ str(ref_iteration) + '/')
							else:
								refTable_d = dar.applyDAR(fitsfile, refTable_tf, plot=True, instrument=osiris, plotdir=plotDir_n + 'dar_a/'+ str(ref_iteration) + '/')
								#print(starlist.info)
								plot_dots(refTable_d,starlist,filename,PA,plotDir_n + 'dots_b/' + str(ref_iteration) + '/')
								# errs = np.hypot(refTable_d['x0e'], refTable_d['x0e'])
								# plot_mags(refTable_d['m0'], errs, filename, plotDir_n + 'mags_b/' + str(ref_iteration) + '/')
			

							try:

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
								
								# trans_list = None

								# print(i, filename, len(refTable_d), 'Reference stars,', len(starlist), 'OSIRIS stars,')

								# starlist_corrected = Table(starlist,copy=True)
								starlist_corrected = starlist[:]

								if current_distortion_correction is not None:
									xt, yt = current_distortion_correction.evaluate(starlist_corrected['x'],starlist_corrected['y'])
									starlist_corrected['x'] = xt 
									starlist_corrected['y'] = yt 




								if single_fit:
									# tab1, tform_inv, tform = gaia_transform(refTable_d, [starlist],trans_list, note='b s-{} {}'.format(i,ref_iteration))
									
									refTable_d_f = refTable_d #.filled(). not filling.
									msc = align.MosaicToRef(
										refTable_d_f, [starlist_corrected], iters=3,
									    dr_tol=rad_tolerance, dm_tol=mag_tolerance_dict[night],
									    outlier_tol=[None, None, None],
										# trans_class=transforms.PolyTransform,
										trans_input=trans_list,
										trans_class=transforms.four_paramNW,
										trans_args=[{'order': 1}, {'order': 1}, {'order': 1}],
									    use_vel=False,
									    use_ref_new=False,
									    update_ref_orig=False, 
									    mag_trans=True,
									    mag_lim=mag_limits, #[6,16],
									    weights='both,std',
									    calc_trans_inverse=True,    
									    init_guess_mode='miracle', verbose=0)

									#trans_class=transforms.PolyTransform,
									#trans_args=[{'order': 1}, {'order': 1}, {'order': 1}], 
									# trans_class=transforms.four_param,
									msc.fit()
									# tab1 = msc.ref_table
									tform = msc.trans_list
									tform_inv = msc.trans_list_inverse

									# if i == 0:
									# 	print('tform px py', tform[0].px, tform[0].py)


									#I only want to save the transformation from the subsequent iterations, not tab1.
									if fp_iteration == 0:
										tab1_initial[night].append(msc.ref_table)
										tab1 = msc.ref_table
									elif fp_iteration > 0:
										# print(tab1_initial)
										# print(tab1_initial[0])
										tab1 = tab1_initial[night][i]
										# tab2 = msc.ref_table


									j = 0
								else:
									j=i  # i=osiris index, j=gaia index


								ref_idx = np.where(tab1['ref_orig'] == True)[j]

								if len(ref_idx) <= 10:
									print(i, filename, 'Only', len(ref_idx),'matched, skipping')
									errorlist.append(filename[:-4])
									continue

								with open(tform_file_1, 'wb') as temp:
									pickle.dump(tform, temp)

								
								# last_good_transform = tform


								print(i, filename, len(refTable_d), 'Reference stars,', len(starlist_corrected), 'OSIRIS stars,', len(ref_idx), 'matches')
								# ids1 = np.where(tab1['name_in_list'] == 's1')[0]
								Ox = tab1['x_orig'][ref_idx,j]
								Oy = tab1['y_orig'][ref_idx,j]
								GRa = tab1['x0'][ref_idx]
								GDec = tab1['y0'][ref_idx]
								# ids = gaiaData['source_id'][ref_idx]   #ID's of gaia sources
								# ids = refTable_d['source_id'][ref_idx]   #ID's of gaia sources
								# ids = refTable_d['name'][ref_idx]   #ID's of stars?   I lose the catalogue names in MosaicToRef, so can't use this.
								ids = tab1['name'][ref_idx]   #ID's of stars?



								px = tform_inv[j].px
								py = tform_inv[j].py
								theta = math.atan2(px[2],px[1])
								scale = math.cos(theta) / px[1]

								used_files.append(filename)
								scales.append(scale)
								rotations.append(math.degrees(theta))
								PAs.append(PA)

								if not recalculate_plate_scale:

									# for k, x in enumerate(Gx):
									# 	Gx[k], Gy[k] = tform_inv[j].evaluate(Gx[k],Gy[k])	#transform Hubble reference into pixels with 4-param transform.
									
					
									Gx, Gy = tform_inv[j].evaluate(GRa,GDec)  #should work fine? Don't need the loop above?
									
									ORa, ODec = tform[j].evaluate(Ox,Oy)


									w = tab1['w'][ref_idx]
									print('Number of times weight=0.0:', (w==0.0).sum())

									# print(np.median(w))
									x_O.extend(Ox)
									y_O.extend(Oy)
									x_G.extend(Gx)
									y_G.extend(Gy)
									weights.extend(w)
									idnums.extend(ids)
									starlist_num.extend([i] * len(ids))
									night_col.extend([night] * len(ids))

									xe_O.extend(tab1['xe_orig'][ref_idx,j])
									ye_O.extend(tab1['ye_orig'][ref_idx,j])
									xe_G.extend(tab1['x0e'][ref_idx])
									ye_G.extend(tab1['y0e'][ref_idx])
									used_in_trans_1.extend(tab1['use_in_trans'][ref_idx])

									Ra_O.extend(ORa)
									Dec_O.extend(ODec)
									Ra_G.extend(GRa)
									Dec_G.extend(GDec)

									plot_image(Gx,Gy,Ox,Oy,fitsfile,plotDir_n + 'img/'+ str(ref_iteration) + '/',tab1['use_in_trans'][ref_idx])
									plot_image_dots(Gx,Gy,Ox,Oy,fitsfile,plotDir_n + 'img_d/'+ str(ref_iteration) + '/',tab1['use_in_trans'][ref_idx])
									plot_quiver(Ox,Oy,Gx,Gy,filename[:-4],plotDir_n + 'quiver/'+ str(ref_iteration) + '/',tab1['use_in_trans'][ref_idx])
									if show_plots:
										plt.show()

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
						# plt.close('all')
						# quit()


				#--------------------------------------------Night loop finishes here----------------------------------------------


				print('--------------------')
				print('Succeeded for', len(successlist), 'files.')
				print('Failed for', len(errorlist), 'files:')
				print(errorlist)

				print('Mean scale =', np.mean(scales))
				offset = np.array(rotations)-np.array(PAs)
				print('Mean rotation offset =', np.mean(offset))

				if recalculate_plate_scale:
					transform_table = Table([used_files, scales, rotations, PAs, offset], names = ('File','Scale', 'Rotation','PA','Difference'))
					ascii.write(transform_table,'{}t_params_{}_{}.txt'.format(resultDir,ref_iteration,night),format='fixed_width', overwrite=True)
				else:
					x_O = np.array(x_O)
					y_O = np.array(y_O)
					x_G = np.array(x_G)
					y_G = np.array(y_G)
					weights = np.array(weights)
					idnums = np.array(idnums)
					starlist_num = np.array(starlist_num)
					night_col = np.array(night_col)
					used_in_trans_1 = np.array(used_in_trans_1)

					output_table = Table([x_O,y_O,x_G,y_G,weights,idnums,xe_O,ye_O,xe_G,ye_G,Ra_O,Dec_O,Ra_G,Dec_G, night_col,starlist_num,used_in_trans_1], names = ('x_OSIRIS','y_OSIRIS','x_REF','y_REF','Weight','REF_ID','xe_OSIRIS','ye_OSIRIS','xe_REF','ye_REF','Ra_OSIRIS','Dec_OSIRIS','Ra_REF','Dec_REF','Night','Frame','UiT'),)
					ascii.write(output_table,dist_file,format='fixed_width', overwrite=True)

					# distortion_data = output_table
					distortion_data = ascii.read(dist_file,format='fixed_width')
					#I generate this for each 4p loop and overwrite it.




			# use the parameters in output table to calculate a new distortion solution (check fit_legendre script)0

			# use that distortion solution to update correction_list

			# used_in_trans_B = distortion_data['UiT'].data
			used_in_trans_B = np.array(distortion_data['UiT'].data) #copy of data, not pointer to
			# print(used_in_trans_B)
			# print(type(used_in_trans_B[0]),used_in_trans_B[0])
			used_in_trans_B = (used_in_trans_B == "True")  #Converts "True" strings to True booleans
			used_in_trans_B.tolist()

			#works fine if loading previous, fails if using new.  wtf. Must be reading the wrong onw? read in to guarantee match.
			# print(used_in_trans_B)
			# print(type(used_in_trans_B[0]),used_in_trans_B[0])
			if trim_not_used_in_trans:
				print('Trimming distortion_list to only include stars used in transformation')
				distortion_data = distortion_data[used_in_trans_B]


			x = distortion_data['x_OSIRIS'].data
			y = distortion_data['y_OSIRIS'].data
			xref = distortion_data['x_REF'].data
			yref = distortion_data['y_REF'].data
			weights = distortion_data['Weight'].data
			gaia_id = distortion_data['REF_ID'].data
			xe = distortion_data['xe_OSIRIS'].data
			ye = distortion_data['ye_OSIRIS'].data
			xeref = distortion_data['xe_REF'].data
			yeref = distortion_data['ye_REF'].data
			frame = distortion_data['Frame'].data
			night_c = distortion_data['Night'].data
			
			ra = distortion_data['Ra_OSIRIS'].data
			dec = distortion_data['Dec_OSIRIS'].data
			raref  = distortion_data['Ra_REF'].data
			decref  = distortion_data['Dec_REF'].data

			# order = 5 #5
			degree = order

			# plot_image(Gx,Gy,Ox,Oy,fitsfile,plotDir_n + 'img/'+ str(ref_iteration) + '/',tab1['use_in_trans'][ref_idx])
			# plot_image_dots(Gx,Gy,Ox,Oy,fitsfile,plotDir_n + 'img_d/'+ str(ref_iteration) + '/',tab1['use_in_trans'][ref_idx])
			# plot_quiver(Ox,Oy,Gx,Gy,filename[:-4],plotDir_n + 'quiver/'+ str(ref_iteration) + '/',tab1['use_in_trans'][ref_idx])
			# if show_plots:
			# 	plt.show()

			# for i, w in enumerate(weights):
			# 	if w == 0.0:
			# 		fix = 1/np.sqrt((xe[i]*0.01)**2 + (ye[i]*0.01)**2 + xeref[i]**2 + yeref[i]**2)
			# 		weights[i] = fix
			# 		print('Correcting weight {} at line {} to be {}. (from outlier rejection)'.format(w,i+2,fix))

			ind_r = np.where(weights == 0.0)[0] #index of stars rejected with outlier_rejection, have their weights set to 0.0
			weights[ind_r] = 1/np.sqrt((xe[ind_r]*0.01)**2 + (ye[ind_r]*0.01)**2 + xeref[ind_r]**2 + yeref[ind_r]**2)
			print('{} stars rejected due to outlier_rejection'.format(len(ind_r)))

			print('Finding outliers')
			outliers, bad_names, m_distances = find_outliers(x,y,xref,yref,gaia_id)
			print('Outliers found')
			# if solution == '2021':
			# 	# manual_outliers = [3536, 3726, 3917, 4687, 4885, 5076]
			# 	# manual_outliers = [3538, 3728, 3919, 4689, 4887, 5078]
			# 	manual_outliers = []
			# 	outliers[manual_outliers] = True

			print('Mean of mahalanobis distances:', np.mean(m_distances))
			print('Standard deviation of mahalanobis distances:', np.std(m_distances))

			bad_names.sort()
			bad_count = Counter(bad_names)

			all_count = Counter(gaia_id)
			bad_stars = [k for (k,v) in bad_count.items() if v > 11]
			# print(bad_count)
			# print('Bad stars with >11 outliers:', bad_stars)
			# plt.bar(*zip(*bad_count.items()))
			# plt.bar(range(len(bad_count)),bad_count.values())
			# plt.show()
			
			if False:
				print('Stars with >10% bad observations (at least 2):')
				for k,v in bad_count.items():
					a = all_count[k]
					if v/a > 0.1 and v > 1:
						print(k,v, 'bad out of', a, '=', v/a)



			include = ~outliers 


			# include = slice(0,None)    #outliers already removed in distortion.py. Can remove all references to include.

			# transModel = transforms.LegClipSplineTransform(x, y, xref, yref, degree, weights=None, kx=None, ky=None, niter=0, sigma=3)
			# tform_leg = transforms.LegTransform.derive_transform(x[include], y[include], xref[include], yref[include], order, m=None, mref=None,init_gx=None, init_gy=None, weights=weights[include], mag_trans=True)
			print('Fitting Legendre Polynomial')
			tform_leg = transforms.LegTransform.derive_transform(x[include], y[include], xref[include], yref[include], order, m=None, mref=None,init_gx=None, init_gy=None, weights=None, mag_trans=True)
			#   Defines a bivariate legendre tranformation from x,y -> xref,yref using Legnedre polynomials as the basis.



			current_distortion_correction = tform_leg



			xts,yts = tform_leg.evaluate(x,y)
			# unexplained_x_error = 0.005
			# unexplained_y_error = 0.005
			unexplained_x_error = 0
			unexplained_y_error = 0
			if centred:
				a,b = tform_leg.evaluate(1024,1024)
				x_central = a-1024
				y_central = b-1024
				xts = xts - x_central
				yts = yts - y_central
				xref = xref - x_central
				yref = yref - y_central

			# plot_quiver(xref[include], yref[include], xts[include], yts[include], order)


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

			# median_4p_residuals.append(np.median(residuals_b[include]))
			# with open(resultDir + 'fp_iteration_residuals.txt', 'w') as temp:
			# 	for r, resid in enumerate(median_4p_residuals):
			# 		# temp.write(str(resid) + '\n')




			x_coefficient_names = []
			x_coefficient_values = []
			y_coefficient_names = []
			y_coefficient_values = []

			for param in tform_leg.px.param_names:
				# print(getattr(tform_leg.px, param))
				a = getattr(tform_leg.px, param)
				x_coefficient_names.append(a.name)
				x_coefficient_values.append(a.value)


			for param in tform_leg.py.param_names:
				a = getattr(tform_leg.py, param)
				y_coefficient_names.append(a.name)
				y_coefficient_values.append(a.value)

			output_table = Table([x_coefficient_names,x_coefficient_values, y_coefficient_names,y_coefficient_values], names = ('px_name','px_val','py_name','py_val'),)
			ascii.write(output_table,'{}distortion_coefficients_{}_{}_{}.txt'.format(resultDir,solution_year,ref_iteration,fp_iteration),format='fixed_width', overwrite=True)




		#------------------------fp_iteration ends here---------------------------
	


		#s


		print('len correction list',len(correction_list), 'ref_iteration', ref_iteration)
		print(correction_list)
		if len(correction_list) == ref_iteration:
			correction_list.append(tform_leg) 
		else:
			print('Correction list should have length', ref_iteration, 'but has length', len(correction_list))
			quit()





		if centred:
			xref = xref + x_central 		# x_ref is already centred, so un-centre it.
			yref = yref + y_central			
			xc, yc = correction_list[ref_iteration].evaluate(1024,1024)
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

		plt.figure(num=1,figsize=(6,6),clear=True)
		ptsize = 0.5
		plt.scatter(x[include],y[include],s=ptsize,edgecolors='none',label='OSIRIS')
		plt.scatter(xref[include], yref[include],s=ptsize,edgecolors='none',label='REF')
		plt.scatter(xts[include],yts[include],s=ptsize,edgecolors='none',label='Legendre')
		# plt.scatter(gx,gy,ptsize,c='r',alpha=0.6,label='Gaia')
		plt.xlim(-400,2448)
		plt.ylim(-400,2448)
		plt.title('Scatter ' + solution_year)
		plt.legend(loc='upper right')
		plt.savefig('{}scatter_{}_{}.pdf'.format(plotDir,solution_year,ref_iteration), bbox_inches='tight')
		# plt.show()

		quiv_scale=200
		# plt.close('all')
		plt.figure(num=1,figsize=(6,6),clear=True)
		q = plt.quiver(x1,y1,(x2-x1),(y2-y1), color='black', scale=quiv_scale, angles='xy',width=0.003)
		quiv_label_val = 10.0
		quiv_label = '{} pix'.format(quiv_label_val)
		plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
		plt.xlim(-400,2448)
		plt.ylim(-400,2448)
		# plt.axis('equal')
		# plt.set_aspect('equal','box')
		plt.title('Order ' + str(order) + ' legendre fit')
		plt.savefig('{}quiver_legendre_{}_{}.pdf'.format(plotDir,solution_year,ref_iteration), bbox_inches='tight',dpi=200)
		plt.savefig('{}quiver_legendre_{}_{}.jpg'.format(plotDir,solution_year,ref_iteration), bbox_inches='tight',dpi=200)


		# quiv_scale=50
		# plt.close('all')
		plt.figure(num=1,figsize=(6,6),clear=True)
		q = plt.quiver(x1,y1,(x2_c-x1),(y2_c-y1), color='black', scale=quiv_scale, angles='xy',width=0.003)
		# quiv_label = '10 pix'
		# quiv_label_val = 10.0
		plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
		plt.xlim(-400,2448)
		plt.ylim(-400,2448)
		# plt.axis('equal')
		# plt.set_aspect('equal','box')
		plt.title('Order ' + str(order) + ' legendre fit (centred)')
		plt.savefig('{}quiver_legendre_c_{}_{}.jpg'.format(plotDir,solution_year,ref_iteration), bbox_inches='tight',dpi=200)

		# quiv_scale=200
		# plt.close('all')
		plt.figure(num=1,figsize=(6,6),clear=True)
		q = plt.quiver(x[include],y[include],(xref[include]-x[include]),(yref[include]-y[include]), color='black', scale=quiv_scale, angles='xy',width=0.0004)
		q = plt.quiver(x[outliers],y[outliers],(xref[outliers]-x[outliers]),(yref[outliers]-y[outliers]), color='red', scale=quiv_scale, angles='xy',width=0.0004)
		# quiv_label = '10 pix'
		# quiv_label_val = 10.0
		plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
		plt.xlim(-400,2448)
		plt.ylim(-400,2448)
		# plt.axis('equal')
		# plt.set_aspect('equal','box')
		plt.title('All matched stars')
		plt.savefig('{}quiver_all_{}_{}.pdf'.format(plotDir,solution_year,ref_iteration), bbox_inches='tight',dpi=200)
		plt.savefig('{}quiver_all_{}_{}.jpg'.format(plotDir,solution_year,ref_iteration), bbox_inches='tight',dpi=200)

		quiv_scale=200
		# plt.close('all')
		plt.figure(num=1,figsize=(6,6),clear=True)
		q = plt.quiver(x[include],y[include],(xref[include]-x_central-x[include]),(yref[include]-y_central-y[include]), color='black', scale=quiv_scale, angles='xy',width=0.0004)
		q = plt.quiver(x[outliers],y[outliers],(xref[outliers]-x_central-x[outliers]),(yref[outliers]-y_central-y[outliers]), color='red', scale=quiv_scale, angles='xy',width=0.0004)
		# quiv_label = '10 pix'
		# quiv_label_val = 10.0
		plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
		plt.xlim(-400,2448)
		plt.ylim(-400,2448)
		# plt.axis('equal')
		# plt.set_aspect('equal','box')
		plt.title('All matched stars (centred)')
		plt.savefig('{}quiver_all_c_{}_{}.pdf'.format(plotDir,solution_year,ref_iteration), bbox_inches='tight',dpi=200)
		plt.savefig('{}quiver_all_c_{}_{}.jpg'.format(plotDir,solution_year,ref_iteration), bbox_inches='tight',dpi=200)

		plt.figure(num=1,figsize=(6,6),clear=True)
		plt.hist(dr,bins=20, linewidth=1,edgecolor='black')
		plt.xlabel('Distortion (pixels)')
		plt.xticks([0,5,10,15,20])
		plt.title('Distortion Histogram')
		plt.savefig('{}dist_histogram_{}_{}.jpg'.format(plotDir,solution_year,ref_iteration), bbox_inches='tight',dpi=200)

		# ------------------- plots from find_legendre_order.py--------------------------------

		

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


		if debugging:
			k = np.where(distortion_data['REF_ID'] == debug_star_id)
			print('k',k)
			k = np.where(gaia_id == debug_star_id)
			print('k',k)

			# print(distortion_data[k])
			# print(distortion_data['REF_ID'][k])
			for l in distortion_data[k]:
				print(l)
			with open(resultDir + 'debug_output.txt', 'a') as temp:
				temp.write(f'Reference b: Star {debug_star_id} found on line {k}\n')
				for m, l in enumerate(distortion_data[k]):
					temp.write(f'{l}\n')
					temp.write(f'x_t:   {xts[k][m]:.3f}  y_t:   {yts[k][m]:.3f}\n')
					temp.write(f'x_ref: {xref[k][m]:.3f}. y_ref: {yref[k][m]:.3f}\n')
					temp.write(f'Residual_b: {residuals_b[k][m]:.6f}\n')
					temp.write('\n')
				# temp.write(f'Residuals: {residuals[k]}\n')
				temp.write('\n\n\n')





		#----------------------------------




		plt.figure(num=3,figsize=(8,6),clear=True)
		plt.hist(residuals_b_x[include],bins=50,histtype='step',range=(-50,50),label='x')
		plt.hist(residuals_b_y[include],bins=50,histtype='step',range=(-50,50),linestyle=('dashed'),label='y')
		# if reference_instrument == 'GAIA':
		# 	plt.ylim([0,1400])
		plt.xlabel('Distortion Residual / standard deviation in measurements')
		plt.ylabel('Count')
		plt.legend()
		plt.title('Distortion residuals_b histogram')
		plt.savefig(plotDir + 'residuals_b_std_' + str(ref_iteration) + '.jpg', bbox_inches='tight',dpi=200)

		plt.figure(num=3,figsize=(8,6),clear=True)
		plt.hist(residuals_b[include],bins=20,range=(0,3))
		# if reference_instrument == 'GAIA':
		# 	plt.ylim([0,1400])
		plt.xlabel('Residual (pixels)')
		plt.ylabel('Count')
		plt.title('Transformation residuals_b histogram')
		plt.savefig(plotDir + 'residuals_b_' + str(ref_iteration) + '.jpg', bbox_inches='tight',dpi=200)

		plt.figure(num=3,figsize=(8,6),clear=True)
		plt.hist(1/(weights*0.01)[include],bins=20,range=(0,0.5))
		# if reference_instrument == 'GAIA':
		# 	plt.ylim([0,1400])
		plt.xlabel('Uncertainty (pixels)')
		plt.ylabel('Count')
		plt.title('Uncertainty histogram')
		plt.savefig(plotDir + 'uncertainties_' + str(ref_iteration) + '.jpg', bbox_inches='tight',dpi=200)

		plt.figure(num=3,figsize=(8,6),clear=True)
		plt.hist(distances_tr[include],bins=20,range=(0,25))
		# if reference_instrument == 'GAIA':
		# 	plt.ylim([0,1400])
		plt.xlabel('Correction (pixels)')
		plt.ylabel('Count')
		plt.title('Correction size histogram')
		plt.savefig(plotDir + 'correction_' + str(ref_iteration) + '.jpg', bbox_inches='tight',dpi=200)

		plt.figure(num=3,figsize=(8,6),clear=True)
		plt.scatter(distances_tr[include],residuals_b[include],4,alpha=0.6,label='Points')
		# plt.scatter(distances_tr[~include],residuals[~include],4,c='r',alpha=0.6,label='Outliers')
		if reference_instrument == 'GAIA':
			# plt.xlim([0,14])
			# plt.ylim([0,13])
			pass
		elif reference_instrument == 'Hubble':
			plt.xlim([0,25])
			plt.ylim([0,5])
			# pass
		plt.xlabel('Correction size (pixels)')
		plt.ylabel('Residual (pixels)')
		plt.title('Residuals_b vs correction')
		plt.savefig(plotDir + 'r_v_c_' + str(ref_iteration) + '.jpg', bbox_inches='tight',dpi=200)		

		plt.figure(num=3,figsize=(8,6),clear=True)
		plt.scatter(distances_tr[include],residuals_b_std[include],4,alpha=0.6,label='Points')
		# plt.scatter(distances_tr[~include],residuals[~include],4,c='r',alpha=0.6,label='Outliers')
		if reference_instrument == 'GAIA':
			# plt.xlim([0,14])
			# plt.ylim([0,13])
			pass
		elif reference_instrument == 'Hubble':
			# plt.xlim([0,25])
			# plt.ylim([0,3])
			pass
		plt.xlabel('Correction size (pixels)')
		plt.ylabel('Residual (standard deviations)')
		plt.title('Residuals_b vs correction')
		plt.savefig(plotDir + 'rvc_std' + str(ref_iteration) + '.jpg', bbox_inches='tight',dpi=200)	

		# plt.figure(num=3,figsize=(8,6),clear=True)
		# plt.scatter(np.hypot(x_std,y_std)[include],residuals[include],4,alpha=0.6,label='Points')
		# # plt.scatter(distances_tr[~include],residuals[~include],4,c='r',alpha=0.6,label='Outliers')
		# if reference_instrument == 'GAIA':
		# 	plt.xlim([0,14])
		# 	plt.ylim([0,13])
		# elif reference_instrument == 'Hubble':
		# 	# plt.xlim([0,25])
		# 	# plt.ylim([0,3])
		# 	pass
		# plt.xlabel('Bootstrap Standard deviation (pixels)')
		# plt.ylabel('Residual (pixels)')
		# plt.title('Std dev vs residuals')
		# plt.savefig(plotDir + 'rvstd' + str(order) + '.jpg', bbox_inches='tight',dpi=200)	


		plt.figure(num=3,figsize=(8,6),clear=True)
		plt.scatter(distances_ref[include],residuals_b[include],4,alpha=0.6,label='Points')
		if reference_instrument == 'Hubble':
			plt.xlim([0,25])
			plt.ylim([0,5])		
		plt.xlabel('Distortion size (pixels)')
		plt.ylabel('Residual (pixels)')
		plt.title('Residuals_b vs Distortions')
		plt.savefig(plotDir + 'r_v_d_' + str(ref_iteration) + '.jpg', bbox_inches='tight',dpi=200)	




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








		# --------------------------------------- Section A ---------------------------------------------------------------

		if len(refTable_current) == 0:
			with open(refTable_current_filename, 'rb') as temp:
				refTable_current = pickle.load(temp) 

		#------------------------Night loop 1 starts here---------------------------
		for night in obs_nights:
			if night == 'n1':
				target = 'm15'
				nightDir = '/u/mfreeman/work/d/n1/'
				cleanDir = nightDir + 'clean/m15_kn3_tdOpen/'
				hubble_file = '/g/lu/data/m15/hst_ref/NGC7078cb.pm' #'je0o61lzq_flt.xymrduvqpk'
				targetID = 1745948323734090368   #gaia ID of target star. Ra and Dec used in prepare_gaia_for_flystar()
				# ra_field = '21:29:57.60'  #approximate centre of FoV, selects gaia stars in radius  
				# dec_field = '12:10:28.3'
				ra_field = 322.48999069
				dec_field = 12.17453385
				radius = 20   #arcseconds
				minmag = 15.4  #dimmest mag for cut
				single_fit = True  #run Mosaic2Ref for each image individually
				rad_tolerance = [0.4, 0.4, 0.2]
				mag_tolerance = [2, 2, 2]
				mag_limits = [6,16]
				bad_files = ['ci200804_a022007_flip_0.8_stf.lis',
							'ci200804_a026012_flip_0.8_stf.lis',
							'ci200804_a027003_flip_0.8_stf.lis',	
							]
				dont_trim = ['ci200804_a014004_flip_0.8_stf.lis', 'ci200804_a026009_flip_0.8_stf.lis']
				# 026002 is bad either way. Didn't find the brightest star. Could be that flag I selected. Try without?
				sl = slice(0,None)
				# sl = slice(34,38)
				show_plots = False

			elif night == 'n2':
				target = 'm92'
				nightDir = '/u/mfreeman/work/d/n2/'
				cleanDir = nightDir + 'clean/m92_kp_tdOpen/'
				hubble_file = '/g/lu/data/m92/hst_ref/NGC6341cp.pm' #'idk901xpq_flt.xymrduvqpk'
				targetID = 1360405503461790848  #gaia ID of target star. Ra and Dec used in prepare_gaia_for_flystar()
				ra_field = 259.285096306648     #approximate centre of FoV, selects gaia stars in radius
				dec_field = 43.13751895071527
				radius = 20   #arcseconds
				minmag = 15.4  #dimmest mag for cut
				rad_tolerance = [0.4, 0.4, 0.2]
				mag_tolerance = [2, 2, 2]
				mag_limits = None
				single_fit = True  #run Mosaic2Ref for each image individually
				bad_files = []
				dont_trim = []
				sl = slice(0,None)
				# sl = slice(151,None)
				show_plots = False    

			elif night == 'n3':
				target = 'm92'
				nightDir = '/u/mfreeman/work/d/n3/'
				cleanDir = nightDir + 'clean/m92_kp_tdOpen/'
				hubble_file = '/g/lu/data/m92/hst_ref/NGC6341cp.pm' #'idk901xpq_flt.xymrduvqpk'
				targetID = 1360405503461790848  #gaia ID of target star. Ra and Dec used in prepare_gaia_for_flystar()
				ra_field = 259.285096306648     #approximate centre of FoV, selects gaia stars in radius
				dec_field = 43.13751895071527
				radius = 20   #arcseconds
				minmag = 15.4  #dimmest mag for cut
				rad_tolerance = [0.4, 0.4, 0.2]
				mag_tolerance = [2, 2, 2]
				mag_limits = None
				single_fit = True  #run Mosaic2Ref for each image individually
				bad_files = []
				dont_trim = []
				sl = slice(0,None)
				# sl = slice(151,None)
				show_plots = False    

			elif night == 'n4':
				target = 'm92'
				nightDir = '/u/mfreeman/work/d/n4/'
				cleanDir = nightDir + 'clean/m92_kp_tdOpen/'
				hubble_file = '/g/lu/data/m92/hst_ref/NGC6341cp.pm' #'idk901xpq_flt.xymrduvqpk'
				targetID = 1360405503461790848  #gaia ID of target star. Ra and Dec used in prepare_gaia_for_flystar()
				ra_field = 259.285096306648     #approximate centre of FoV, selects gaia stars in radius
				dec_field = 43.13751895071527
				radius = 20   #arcseconds
				minmag = 15.4  #dimmest mag for cut
				rad_tolerance = [0.4, 0.4, 0.2]
				mag_tolerance = [2, 2, 2]
				mag_limits = None
				single_fit = True  #run Mosaic2Ref for each image individually
				bad_files = ['ci200814_a032008_flip_0.8_stf.lis']
				dont_trim = []
				sl = slice(0,None)
				# sl = slice(151,None)
				show_plots = False  

			elif night == 'n5':
				target = 'm15'
				nightDir = '/u/mfreeman/work/d/n5/'
				cleanDir = nightDir + 'clean/m15_kn3_tdhBand/'
				hubble_file = '/g/lu/data/m15/hst_ref/NGC7078cb.pm' #'icbe05m9q_flt.xymrduvqpk'     #there are other files too?
				targetID = 1745948328028761984  #gaia ID of target star. Ra and Dec used in prepare_gaia_for_flystar()
				ra_field = 322.4912419147502  #approximate centre of FoV, selects gaia stars in radius
				dec_field = 12.164721331771977	
				radius = 20   #arcseconds
				minmag = 15.4  #dimmest starlist mag for cut
				rad_tolerance = [0.4, 0.4, 0.2]
				mag_tolerance = [2, 2, 2]
				mag_limits = [6,16]
				single_fit = True  #run Mosaic2Ref for each image individually
				bad_files = ['ci211024_a011009_flip_0.8_stf.lis','ci211024_a012010_flip_0.8_stf.lis','ci211024_a020012_flip_0.8_stf.lis']
				dont_trim = []
				sl = slice(0,None)
				# sl = slice(79,None)
				show_plots = False 
			else:
				print('No night selected')
				quit()

			starfindDir = cleanDir + 'starfinder/'
			stackDir = cleanDir + 'stacks/'
			plotDir_n = plotDir + night + '/'

			osiris_filenames = get_osiris_files(stackDir)

			sl2 = slice(0,None) 
			osiris_filenames = osiris_filenames[sl2]
			print(osiris_filenames)
			print(len(osiris_filenames), 'OSIRIS images')




			# combined_ref_filename = resultDir + 'combined_ref_table_' + str(ref_iteration) + '.txt'
			combined_ref_filename = '{}combined_ref_table_{}_{}.txt'.format(resultDir,night,ref_iteration)
			combined_ref_filename_previous = previous_results_location + reference_instrument + '_' + fitmode + '/' + 'combined_ref_table_' + str(ref_iteration) + '.txt'
			combined_ref_filename_previous = '{}{}_{}/combined_ref_table_{}_{}.txt'.format(previous_results_location,reference_instrument,fitmode,night,ref_iteration)

			if os.path.exists(combined_ref_filename):
				combined_ref_filename_toUse = combined_ref_filename
			elif os.path.exists(combined_ref_filename_previous):
				combined_ref_filename_toUse = combined_ref_filename_previous
			else:
				combined_ref_filename_toUse = None

			# if (not create_combined_reflist) and os.path.exists(combined_ref_filename):
			# 	print('a: Loading previous combined refTable')
			# 	with open(combined_ref_filename, 'rb') as combined_ref_file:
			# 		refTable_a = pickle.load(combined_ref_file)
			
			if (not create_combined_reflist) and (combined_ref_filename_toUse is not None):
				print('a: Loading previous combined refTable')
				with open(combined_ref_filename_toUse, 'rb') as combined_ref_file:
					refTable_a = pickle.load(combined_ref_file)

			else:
				list_of_starlists = []
				print('a: Generating new combined refTable')

				for i, filename in enumerate(osiris_filenames):
					if filename in bad_files:
						print('{} {} flagged as bad, skipping'.format(i,filename))
						# errorlist.append(filename[:-4])
					else:
						print('{} {} applying distortion correction'.format(i,filename))
						starlist = load_osiris_file(stackDir ,filename)
						# plt.close('all')

						# ido = np.where(starlist['m'] < 15.5)
						# starlist = starlist[ido]
						fitsfile = cleanDir + filename[:-12] + '.fits'

						PA = get_PA(fitsfile)
						# print('PA', PA)
						# starlist = brightest_n(starlist,170)
						starlist = mag_cut(starlist,0,minmag)
						if not filename in dont_trim:
							starlist = edge_cut(starlist,5)

						if len(starlist) == 0:
							print(i,filename, '0 stars remaining after edge cut, skipping image')
							errorlist.append(filename[:-4])
							continue

						#----------------------------------------
						# plt.figure()
						# plt.hist(starlist['m'])
						# # plt.scatter(starlist['m'],starlist['vxe'],ptsize,alpha=0.2,label='Ref')
						# plt.show()
						# quit()
						#--------------------------------------

						#apply distortion correction
						xt, yt = correction_list[ref_iteration].evaluate(starlist['x'],starlist['y'])

						# if centred:
						# 	xc, yc = correction_list[ref_iteration].evaluate(1024,1024)
						# 	xc -= 1024
						# 	yc -= 1024
						# else:
						# 	xc = 0
						# 	yc = 0

						#I think I always want the centred version. The DAR is calculated based on the header elevation, so I want the central pixel to be undistorted
						#maybe I just always want it centred. Because that is the actual distortion. And it doesn't effect the model.
						xc, yc = correction_list[ref_iteration].evaluate(1024,1024)
						xc -= 1024
						yc -= 1024						
						starlist['x'] = xt - xc
						starlist['y'] = yt - yc



						# if reference_instrument == 'Hubble':
						# 	refTable_t = refTable
						# elif reference_instrument == 'GAIA':
						# 	refTable_t = trim_gaia(refTable,filename,PA)

						plt.figure(num=4,figsize=(6,6),clear=True)
						if recalculate_plate_scale:
							# refTable_d = dar.applyDAR(fitsfile, refTable_t, plot=False, instrument=osiris, plotdir=plotDir + 'dar/')
							starlist = dar.removeDAR(fitsfile,starlist, plot=False, instrument=osiris, plotdir=plotDir_n + 'dar_r/'+ str(ref_iteration) + '/')

						else:
							# refTable_d = dar.applyDAR(fitsfile, refTable_t, plot=True, instrument=osiris, plotdir=plotDir + 'dar/')
							starlist = dar.removeDAR(fitsfile,starlist, plot=True, instrument=osiris, plotdir=plotDir_n + 'dar_r/'+ str(ref_iteration) + '/')

							#print(starlist.info)
							# plot_dots(refTable_current[night],starlist,filename,PA,plotDir_n + 'dots_a/' + str(ref_iteration) + '/')
							#refTable isn't declared if I load section B. So just don't plot it? I already run plot_dots in section B
						


						list_of_starlists.append(starlist)

						# successlist.append(filename[:-4])

					# plt.close('all')
					# quit()

				print('--------------------')
				print('Completed corrections')
				# print('Succeeded for', len(successlist), 'files.')
				# print('Failed for', len(errorlist), 'files:')
				# print(errorlist)


				print(len(list_of_starlists))
				print([len(j) for j in list_of_starlists])

				# do mosaicSelfref on list_of_starlists to generate a master list.


					

				if reference_instrument == 'Hubble':
					# os.makedirs('./transform_files/hubble_{}'.format(fitmode),exist_ok=True)
					# tform_file_ref = './transform_files/hubble_{}/tform_{}_{}.p'.format(fitmode,ref_iteration,'ref')
					# tform_file_ref_last = './transform_files/hubble_{}/tform_{}_{}.p'.format(fitmode,ref_iteration-1,'ref')
					# tform_file_ref_previous = ['/u/mfreeman/work/d/transform_files/hubble/tform_{}.p'.format(i) for i in osiris_filenames]

					os.makedirs('{}combining_ref/'.format(tformDir),exist_ok=True)
					tform_file_ref = '{}combining_ref/tform_{}_{}.p'.format(tformDir,ref_iteration,night)
					tform_file_ref_2 = '{}combining_ref/tform_{}_{}.p'.format(transform_files_location,fitmode,ref_iteration,night)
					tform_file_ref_last = '{}combining_ref/tform_{}_{}.p'.format(tformDir,fitmode,ref_iteration-1,night)
					tform_file_ref_previous = ['/u/mfreeman/work/d/transform_files/hubble/tform_{}.p'.format(i) for i in osiris_filenames]



				elif reference_instrument == 'GAIA':
					#haven't updated to use new location
					tform_file_ref = './transform_files/gaia_{}/tform_{}_{}.p'.format(fitmode,ref_iteration,'ref')
					tform_file_ref_previous = ['/u/mfreeman/work/d/transform_files/gaia/tform_{}.p'.format(i) for i in osiris_filenames]

				
				use_individual_trans_guesses = False

				# if ref_iteration == 2:


				if os.path.exists(tform_file_ref):
					print('Loading transform from {}'.format(tform_file_ref))
					with open(tform_file_ref, 'rb') as trans_file:
						trans_list_temp = pickle.load(trans_file)

					if len(trans_list_temp) == len(list_of_starlists):
						print('Transform guess has correct number of lists {}'.format(len(trans_list_temp)))
						trans_list = trans_list_temp
					else:
						use_individual_trans_guesses = True
		
				elif os.path.exists(tform_file_ref_2):
					print('Loading a previous transform from {}'.format(tform_file_ref_2))
					with open(tform_file_ref_2, 'rb') as trans_file:
						trans_list_temp = pickle.load(trans_file)

					if len(trans_list_temp) == len(list_of_starlists):
						print('Transform guess has correct number of lists: {}'.format(len(trans_list_temp)))
						trans_list = trans_list_temp
					else:
						use_individual_trans_guesses = True

				elif os.path.exists(tform_file_ref_last):
					print('Loading last transform from {}'.format(tform_file_ref_last))
					with open(tform_file_ref_last, 'rb') as trans_file:
						trans_list_temp = pickle.load(trans_file)

					if len(trans_list_temp) == len(list_of_starlists):
						print('Transform guess has correct number of lists: {}'.format(len(trans_list_temp)))
						trans_list = trans_list_temp
					else:
						use_individual_trans_guesses = True

				else:
					use_individual_trans_guesses = True

				if use_individual_trans_guesses:
					print('Loading old transform')
					if os.path.exists(tform_file_ref_previous[0]):
						trans_list = []
						for i, tform_filename in enumerate(tform_file_ref_previous):
							with open(tform_filename, 'rb') as trans_file:
								trans_list.extend(pickle.load(trans_file))					#initial guess for the transformation
																							#each file is a list with one element, so using .extend
					else:
						# trans_list = last_good_transform
						print('No trans_list found')
						trans_list = None



				if fitmode == 'FF':
					set_use_ref_new = False
					set_update_ref_orig = False
				elif fitmode == 'TF':
					set_use_ref_new = True
					set_update_ref_orig = False
				elif fitmode == 'FT':				
					set_use_ref_new = False
					set_update_ref_orig = True
				elif fitmode == 'TT':
					set_use_ref_new = True
					set_update_ref_orig = True	
				else:
					print('Fitmode not recognised: {}'.format(fitmode))

				print('Creating combined reference frame with MosaicToRef')
				msc = align.MosaicToRef(refTable_current[night],
					list_of_starlists, iters=3,
					dr_tol=[0.5,0.5,0.5], dm_tol=[2,2,2],
					outlier_tol=[None,None,None],
					# trans_class=transforms.PolyTransform,
					trans_input=trans_list,
					trans_class=transforms.four_paramNW,
					trans_args=[{'order': 1}, {'order': 1}, {'order': 1}],
					use_vel=use_flystar_velocity,
					mag_trans=True,
					mag_lim=[6,13], #[6,16],
					weights=None,
					use_ref_new = set_use_ref_new,
					update_ref_orig = set_update_ref_orig,
					calc_trans_inverse=True,    
					init_guess_mode='miracle', verbose=0)




				msc.fit()
				# tab2 = msc.ref_table

				# if reference_instrument == 'Hubble':
				# 	tform_file_ref = './transform_files/hubble/tform_{}_{}.p'.format(ref_iteration,'ref')
				# elif reference_instrument == 'GAIA':
				# 	tform_file_ref = './transform_files/gaia/tform_{}_{}.p'.format(ref_iteration,'ref')
				# else:
				# 	print('wrong reference instrument', reference_instrument)
				# 	tform_file_ref = None

				with open(tform_file_ref, 'wb') as temp:
					pickle.dump(msc.trans_list, temp)
				
				print('Completed reference table matching')


				refTable_a = msc.ref_table

				with open(combined_ref_filename, 'wb') as temp:
					pickle.dump(msc.ref_table, temp)


			#----------- end creation of combined_ref_file -------------------


			



			#average the star points manually here? I think they are all saved in columns of refTable_a
			


			refTable_b = Table([refTable_a['name']], names=['name'], masked=False)

			if manually_average_star_positions:
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


			# print(refTable_b['use_in_trans'])
			# print(list(refTable_b['use_in_trans']).count(True))
			# print(list(refTable_b['use_in_trans']).count(False))
			# print(sum(refTable_b['use_in_trans']))
			# print(len(refTable_b['use_in_trans']))

			if trim_not_used_in_trans:
				print('Trimming refTable_b to only include stars used in transformation')
				refTable_b = refTable_b[refTable_a['use_in_trans']]			#using refTable_a due to the problem above with refTable_b

			# print(refTable_b['use_in_trans'])
			# print(list(refTable_b['use_in_trans']).count(True))
			# print(list(refTable_b['use_in_trans']).count(False))
			# print(sum(refTable_b['use_in_trans']))
			# print(len(refTable_b['use_in_trans']))
			# quit()

			bla, idx1, idx2 = np.intersect1d(refTable_b['name'].astype('str'),refTable_current[night]['name'],return_indices=True)


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


			residuals_a = np.hypot(refTable_current[night]['x0'][idx2]-refTable_b['x0'][idx1],refTable_current[night]['y0'][idx2]-refTable_b['y0'][idx1])
			mean_residuals_a.append(np.mean(residuals_a))

			if debugging:
				debug_star_id = '5'
				k = np.where(refTable_b['name'] == debug_star_id)
				print('k',k)
				# print(refTable_b[k])
				# print(refTable_b['name'][k])
				for l in refTable_b[k]:
					print(l)
				with open(resultDir + 'debug_output.txt', 'a') as temp:
					temp.write('------------------------------------Ref Iteration {}------------------------------------\n\n'.format(ref_iteration))
					temp.write(f'Reference a: Star {debug_star_id} found on line {k}\n')
					for l in refTable_b[k]:
						temp.write(f'{l}\n\n')

			#applyDAR requires the refTable to be in arcseconds. But it does just use the relative position. So if I convert everything using the plate scale, it should be fine
			# refTable_t['x0'] = refTable_t['x0'] * osiris_scale
			# refTable_t['y0'] = refTable_t['y0'] * osiris_scale
			#this might break lots of things?
			#plot dots may not work.
			

			# do the same thing as normal. I have the starlists, run the macthing, save the matched stars. Stack together, do distortion. Repeat in loop somehow.



			# print('NOT UPDATING REF TABLE')
			# refTable_current = Table(refTable_b,copy=True)
			refTable_current[night] = Table(refTable_b,copy=True)


		#----------------------------Night loop 2 ends here-----------------------




	#---------------------------------ref_iteration loop ends here------------------------------------


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












# def make_dir(plotDir,ref_iteration=0):
# 	if not os.path.exists(plotDir):
# 		os.makedirs(plotDir)
# 	if not os.path.exists(resultDir):
# 		os.makedirs(resultDir)
# 	if not os.path.exists(stackDir):
# 		os.makedirs(stackDir)	
# 	if not os.path.exists(plotDir_n + 'dar_r/'+ str(ref_iteration) + '/'):
# 		os.makedirs(plotDir_n + 'dar_r/'+ str(ref_iteration) + '/')
# 	if not os.path.exists(plotDir_n + 'dar_a/'+ str(ref_iteration) + '/'):
# 		os.makedirs(plotDir_n + 'dar_a/'+ str(ref_iteration) + '/')

def make_dir(plotDir,ref_iteration=0):
	os.makedirs(plotDir, exist_ok=True)
	os.makedirs(resultDir, exist_ok=True)
	os.makedirs(stackDir, exist_ok=True)
	os.makedirs(tformDir, exist_ok=True)
	os.makedirs(plotDir_n + 'dar_r/'+ str(ref_iteration) + '/', exist_ok=True)
	os.makedirs(plotDir_n + 'dar_a/'+ str(ref_iteration) + '/', exist_ok=True)





def fetch_gaia(ra, dec, radius, targetID, night, use_prev_gaia=False):
	gaiaName = 'gaia3_' + night + '.p'
	ra = Angle(ra * u.deg).to_string(unit=u.hour, decimal=True)
	if use_prev_gaia:
		if os.path.exists(gaiaName):
			with open(gaiaName, 'rb') as temp:
				gaiaData = pickle.load(temp)
		else:
			gaiaData = analysis.query_gaia(ra, dec, radius, 'gaiaedr3')
			with open(gaiaName,'wb') as temp:
				pickle.dump(gaiaData,temp)
	else:
		gaiaData = analysis.query_gaia(ra, dec, radius, 'gaiaedr3')

	print(gaiaData.info)
	# idx = np.where(gaia['pmdec'].mask == True)[0]
	# idx=gaiaData['pmra'].mask
	# gaiaData['pmra'][idx] = 0.0
	# gaiaData['pmdec'][idx] = 0.0
	# gaiaData['pmra_error'][idx] = 0.0
	# gaiaData['pmdec_error'][idx] = 0.0


	# gaiaData = gaiaData.filled()      #is this ok?
	# print(gaiaData['ra'])
	print(len(gaiaData), 'Gaia stars')
	return gaiaData


def fetch_hubble(filename):				#loads hubble data from hst1pass
	# filelist = os.listdir(nightDir)
	# hst1pass_file = fnmatch.filter(filelist,'i*.xymrduvqpk')[0]
	print('Loading', filename)
	# hubble_table = ascii.read(nightDir + hst1pass_file,format='commented_header',guess=False,delimiter='\s',comment='#',header_start=166,data_start=168,) #lines start at 0
	# hubble_table = ascii.read(nightDir + hst1pass_file,format='basic',names=['x','y','m','r','d','u','v','q','p','k'])
	column_names = ['r', 'x_0', 'y_0', 'pmx', 'pmy', '1sig_pmx', '1sig_pmy', 'x_M', 'y_M', 'Delta_x', 'Delta_y', 
					'1sig_pmx_mod', '1sig_pmy_mod', 'm_F606W', 'm_F814W', 'rms_F606W', 'rms_F814W', 'qfit_F606W', 'qfit_F814W', 
					'red_chi2_pmx', 'red_chi2_pmy', '1sig_intrcp_pmx', '1sig_intrcp_pmy', 'time_bsln', '1sig_intrcp_pmx_mod', '1sig_intrcp_pmy_mod', 
					'pm_ref', 'N_found', 'N_used', 'ID', 'delta_pmx', 'delta_pmy']

	hubble_table = ascii.read(filename,format='basic',names=column_names)
	return hubble_table

def prepare_hubble_for_flystar(hubble, ra, dec, target, targets_dict=None, match_dr_max=0.2):   #copied from prepare_gaia_for_flystar(). #year=2006 by default?
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
		ra_diff = (ra_centre-ra)*3600     #degrees to arcseconds
		dec_diff = (dec_centre-dec)*3600   

		manual_ra =  -(0) * osiris_scale    #manual offset to align images.
		manual_dec = -(0) * osiris_scale
		ra_diff += manual_ra
		dec_diff += manual_dec
	else:
		print('{} can only be be m15 or m92'.format(target))

	# cos_dec = np.cos(np.radians(dec))
	# x = (hubble['r'] - ra) * cos_dec * 3600.0   # arcsec
	# y = (hubble['d'] - dec) * 3600.0           # arcsec
	#xe = hubble['ra_error'] * cos_dec / 1e3      # arcsec
	#ye = hubble['dec_error'] / 1e3               # arcsec

	# need to make a column of names. Maybe source ID too? Or can I skip it.

	# hubble['source_id'] = range(len(hubble['x_0']))

	# hubble_new = table.Table([hubble['source_id'].data.astype('S19')], names=['name'], masked=False)
	hubble_new = Table([hubble['ID']], names=['name'], masked=False)

	# hubble_new['x0'] = x * -1.0
	# hubble_new['y0'] = y
	#hubble_new['x0e'] = xe
	#hubble_new['y0e'] = ye
	hubble_new['x0'] = (hubble['x_0']+ra_diff) *-1
	hubble_new['y0'] = hubble['y_0']+dec_diff
	hubble_new['x0e'] = hubble['1sig_intrcp_pmx'] * 40/1000   #40mas/pix to arcseconds
	hubble_new['y0e'] = hubble['1sig_intrcp_pmy'] * 40/1000 

	# Also convert the velocities. Note that gaia PM are already * cos(dec)
	#hubble_new['vx'] = hubble['pmra'].data * -1.0 / 1e3 # asec/yr
	#hubble_new['vy'] = hubble['pmdec'].data / 1e3
	#hubble_new['vxe'] = hubble['pmra_error'].data / 1e3
	#hubble_new['vye'] = hubble['pmdec_error'].data / 1e3
	hubble_new['vx'] = hubble['pmx'] / 1000    #pmx is in mas per year. Converting to arcsec per year
	hubble_new['vy'] = hubble['pmy'] / 1000
	hubble_new['vxe'] = hubble['1sig_pmx'] / 1000     #arcsec per year
	hubble_new['vye'] = hubble['1sig_pmy'] / 1000   

	hubble_new['t0'] = hubble['time_bsln']
	# hubble_new['source_id'] = hubble['ID']

	# Find sources without velocities and fix them up.
	# idx = np.where(hubble['pmy'].mask == True)[0]
	# hubble_new['vx'][idx] = 0.0
	# hubble_new['vy'][idx] = 0.0
	# hubble_new['vxe'][idx] = 0.0
	# hubble_new['vye'][idx] = 0.0

	hubble_new['m0'] = hubble['m_F814W']
	# hubble_new['me'] = 1.09/hubble['phot_g_mean_flux_over_error']
	# hubble_new['parallax'] = hubble['parallax']
	# hubble_new['parallax_error'] = hubble['parallax_error']

	# Set the velocities (and uncertainties) to zero if they aren't measured.
	idx = np.where(np.isnan(hubble_new['vx']) == True)[0]
	hubble_new['vx'][idx] = 0.0
	hubble_new['vxe'][idx] = 0.0
	hubble_new['vy'][idx] = 0.0
	hubble_new['vye'][idx] = 0.0

	hubble_new = hubble_new.filled()  #convert masked colunms to regular columns

	if targets_dict != None:
		for targ_name, targ_coo in targets_dict.items():
			dx = hubble_new['x0'] - (targ_coo[0] * -1.0)
			dy = hubble_new['y0'] - targ_coo[1]
			dr = np.hypot(dx, dy)

			idx = dr.argmin()

			if dr[idx] < match_dr_max:
				hubble_new['name'][idx] = targ_name
				print('Found match for: ', targ_name)

	return hubble_new


def get_starfindStack(starfindDir):
	filelist = os.listdir(starfindDir)
	starfindFiles = fnmatch.filter(filelist,'ci*.lis')
	starfindFiles.sort()
	# starfindStack = {}
	starfindStack = []
	for i,filename in enumerate(starfindFiles):
		with open (starfindDir + filename) as starfindFile:
			starfindTable = starfindFile.read().splitlines()
		# starfindStack[filename] = starfindTable
		# starfindStack[i] = starfindTable
		starfindStack.append(starfindTable)
	return starfindStack, starfindFiles

def get_osiris_files(starfindDir):
	filelist = os.listdir(starfindDir)
	starfindFiles = fnmatch.filter(filelist,'ci*.lis')
	starfindFiles.sort()
	return starfindFiles

def load_osiris_file(starfindDir,filename):
	with open (starfindDir + filename) as starfindFile:
			starfindTable = starfindFile.read().splitlines()
	osirisStarList = starlists.StarList.from_lis_file(starfindTable, error=False)
	return osirisStarList


# def project_pos(refTable, starlists):
# 	#Assume that all osiris images have the same epoch.
# 	osirisT = starlists[0]['t'][0]
# 	dt = osirisT - refTable['t0']
# 	refTable['x0'] = refTable['x0'] + refTable['vx'] * dt
# 	refTable['y0'] = refTable['y0'] + refTable['vy'] * dt
# 	return refTable
#hang on, the transformation is different.
#mabye find_transform_new()?


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

def gaia_transform(refTable, osirisStarLists,trans_list=None):
	refTable_f = refTable #.filled()
	msc = align.MosaicToRef(
		refTable_f, osirisStarLists, iters=3,
	    dr_tol=rad_tolerance, dm_tol=mag_tolerance,
	    outlier_tol=[None, None, None],
		# trans_class=transforms.PolyTransform,
		trans_input=trans_list,
		trans_class=transforms.four_paramNW,
		trans_args=[{'order': 1}, {'order': 1}, {'order': 1}],
	    use_vel=False,
	    use_ref_new=False,
	    update_ref_orig=False, 
	    mag_trans=True,
	    mag_lim=mag_limits, #[6,16],
	    weights='both,std',
	    calc_trans_inverse=True,    
	    init_guess_mode='miracle', verbose=0)

	#trans_class=transforms.PolyTransform,
	#trans_args=[{'order': 1}, {'order': 1}, {'order': 1}], 
	# trans_class=transforms.four_param,
	msc.fit()
	tab1 = msc.ref_table
	tform = msc.trans_list
	tform_inv = msc.trans_list_inverse
	return msc.ref_table, tform_inv, msc.trans_list

def plot_quiver(Ox,Oy,Gx,Gy,filename,directory,used):
	os.makedirs(directory, exist_ok=True)		
	quiv_scale=100
	# plt.close('all')
	plt.figure(num=4,figsize=(6,6),clear=True)
	q = plt.quiver(Ox[~used],Oy[~used],(Gx[~used]-Ox[~used]),(Gy[~used]-Oy[~used]), color='b', scale=quiv_scale, angles='xy',width=0.004,label='Not UiT')
	q = plt.quiver(Ox[used],Oy[used],(Gx[used]-Ox[used]),(Gy[used]-Oy[used]), color='black', scale=quiv_scale, angles='xy',width=0.004,label='UiT')
	quiv_label = '10 pix'
	quiv_label_val = 10.0
	plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
	plt.xlim(-400,2448)
	plt.ylim(-400,2448)
	# plt.axis('equal')
	# plt.set_aspect('equal','box')
	plt.title('Matched stars OSIRIS -> Hubble')
	plt.savefig(directory + filename + '.jpg', bbox_inches='tight')
	# plt.show()
	# plt.close()

# @profile
def plot_image(Gx,Gy,Ox,Oy,fitsfile,directory,used):
	os.makedirs(directory, exist_ok=True)		
	# fitsfile = cleanDir + filename[:-12] + '.fits'
	img = fits.getdata(fitsfile)
	img[np.where(img<1)] = 1
	# plt.close('all')
	plt.figure(num=4,figsize=(6,6),clear=True)
	vmin = 10
	vmax = np.max(img) #*0.75
	norm = LogNorm(vmin, vmax)
	plt.imshow(img, cmap='Greys_r', norm=norm, origin = "lower", )
	# ptsize = 5
	# sc = plt.scatter(Ox[~used],Oy[~used],ptsize,c='b',label='Osiris not UiT')
	# sc = plt.scatter(Ox[used],Oy[used],ptsize,c='g',label='Osiris UiT')
	# plt.scatter(Gx,Gy,ptsize,c='r',alpha=0.6,label='Ref')
	plt.xlim(-400,2448)
	plt.ylim(-400,2448)
	plt.title('OSIRIS Frame')
	# plt.legend(loc='upper right')
	plt.savefig(directory + fitsfile[-26:-5] + '.pdf', bbox_inches='tight')
	# plt.close()
	# plt.show()

def plot_image_dots(Gx,Gy,Ox,Oy,fitsfile,directory,used):
	os.makedirs(directory, exist_ok=True)		
	# fitsfile = cleanDir + filename[:-12] + '.fits'
	img = fits.getdata(fitsfile)
	img[np.where(img<1)] = 1
	# plt.close('all')
	plt.figure(num=4,figsize=(6,6),clear=True)
	vmin = 10
	vmax = np.max(img) #*0.75
	norm = LogNorm(vmin, vmax)
	plt.imshow(img, cmap='Greys_r', norm=norm, origin = "lower", )
	ptsize = 5
	sc = plt.scatter(Ox[~used],Oy[~used],ptsize,c='b',label='Osiris not UiT')
	sc = plt.scatter(Ox[used],Oy[used],ptsize,c='g',label='Osiris UiT')
	plt.scatter(Gx,Gy,ptsize,c='r',alpha=0.6,label='Ref')
	plt.xlim(-400,2448)
	plt.ylim(-400,2448)
	plt.title('Star positions (matched)')
	plt.legend(loc='upper right')
	plt.savefig(directory + fitsfile[-26:-5] + '.pdf', bbox_inches='tight')
	# plt.close()
	# plt.show()



def plot_mags(mag, error,filename,directory):
	os.makedirs(directory, exist_ok=True)		
	ptsize = 5
	plt.figure(num=1,figsize=(6,6),clear=True)
	sc=plt.scatter(mag,error,ptsize,c='b')
	plt.xlabel('Magnitude')
	plt.ylabel('Position error (?)')
	plt.savefig(directory + filename + '.jpg', bbox_inches='tight')
	# plt.close()


def plot_dots(refTable, starlist,filename,PA,directory):
	os.makedirs(directory, exist_ok=True)
	PA = math.radians(PA) - 0.02
	with open (cleanDir + filename[:-12] + '.coo') as coofile:
		cooxy = np.genfromtxt(coofile)
	# print(refTable.info)
	#gaia table already in arcseconds relative to coo star
	# idr = np.nonzero(refTable['m0'] < 16)
	idr = np.nonzero(refTable['m0'] < 30)

	gx = refTable['x0'][idr] / osiris_scale + cooxy[0]
	gy = refTable['y0'][idr] / osiris_scale + cooxy[1]

	for i, x in enumerate(gx):
		gx[i], gy[i] = rotate((cooxy[0],cooxy[1]), (gx[i],gy[i]), -PA)
 
	# ids = np.where(starlist['m'] < 16)
	ids = np.where(starlist['m'] < 30)
	# idg = np.where(refTable['m'] < 30)
	# print(starlist.info)

	# print('used', len(starlist['x']), 'osiris stars')
	fitsfile = cleanDir + filename[:-12] + '.fits'
	img = fits.getdata(fitsfile)
	img[np.where(img<1)] = 1
	# ox = np.array(starlist['x'])
	# oy = np.array(starlist['y'])
	# print(starlist['x'])
	# print(ox)
	#convert gaia to pixels
	# match with  star
	# plt.close('all')
	plt.figure(num=4,figsize=(6,6),clear=True)
	vmin = 10
	vmax = np.max(img) #*0.75
	norm = LogNorm(vmin, vmax)
	# norm = Normalize(vmin, vmax)
	# plt.imshow(data, cmap='Greys', origin='lower', norm=norm, interpolation='nearest')
	plt.imshow(img, cmap='Greys_r', norm=norm, origin = "lower", )
	ptsize = 15
	sc = plt.scatter(starlist['x'][ids],starlist['y'][ids],ptsize,linewidths=1,marker='+',c='b',label='Osiris')
	plt.scatter(gx,gy,ptsize,linewidths=1,marker='x',c='r',alpha=0.6,label='Ref')
	plt.xlim(-400,2448)
	plt.ylim(-400,2448)
	plt.title('Star positions (before matching)')
	plt.legend(loc='upper right')
	plt.savefig(directory + filename + '.pdf', bbox_inches='tight')
	# plt.show()
	# plt.close()



def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy

def brightest_n(starlist,n):
	ind = np.argsort(starlist['m'])[:n]
	return starlist[ind]

def mag_cut(starlist,min,max):
	idx = np.where((starlist['m'] >= min) & (starlist['m'] <= max))
	return starlist[idx]


def trim_gaia(refTable,filename,PA):
	#select gaia points that are inside the frame
	with open (cleanDir + filename[:-12] + '.coo') as coofile:
		cooxy = np.genfromtxt(coofile)
	# print('PA = ', PA)
	if abs(PA-0)<5:
		#S eclect a square, with x and y limits
		xmin=(-100 - cooxy[0]) * osiris_scale
		xmax=(2148 - cooxy[0]) * osiris_scale
		ymin=(-100 - cooxy[1]) * osiris_scale
		ymax=(2148 - cooxy[1]) * osiris_scale
		ind = (refTable['x0'] >= xmin) & (refTable['x0'] <= xmax) & (refTable['y0'] >= ymin) & (refTable['y0'] <= ymax)
	elif abs(PA-90)<5:
		xmax=-(-100 - cooxy[1]) * osiris_scale
		xmin=-(2148 - cooxy[1]) * osiris_scale
		ymin=(-100 - cooxy[0]) * osiris_scale
		ymax=(2148 - cooxy[0]) * osiris_scale
		ind = (refTable['x0'] >= xmin) & (refTable['x0'] <= xmax) & (refTable['y0'] >= ymin) & (refTable['y0'] <= ymax)
	else:
		# Select a circle centred on the frame
		x_centre = (1024 - cooxy[0])*osiris_scale
		y_centre = (1024 - cooxy[1])*osiris_scale
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

def get_PA(filename):
	hdr = fits.getheader(filename)
	# return math.radians(hdr['PA_IMAG'])
	return hdr['PA_IMAG']

def edge_cut(starlist,c):
	x1 = c
	x2 = 2048 - c 
	y1 = c 
	y2 = 2048 - c 
	idx = np.where((starlist['x'] >= x1) & (starlist['x'] <= x2) & (starlist['y'] >= y1) & (starlist['y'] <= y2))
	# print('Cutting', len(starlist)-len(starlist[idx]), 'OSIRIS stars from the edge')
	return starlist[idx]


# def find_outliers(x,y,xref,yref,gaia_id):
# 	#scan grid over field. Flag any points as outliers.  Service(2016) did 205x205 pixel bins. removing 3 sigma outliers.

# 	#1 sort data into bins
# 	#2 in each bin remove outliers.

# 	dx = x - xref
# 	dy = y - yref

# 	# bins = np.array([0,256,256*2,256*3,256*4,256*5,256*6,256*7,256*8,])
# 	# bins = np.array([0,2048])
# 	# bins = np.array([0,512,1024,1536,2048])   #looks like 512x512 bins are the best.
# 	bins = np.arange(0,2048+1,256)



# 	xbin = np.digitize(x,bins)
# 	ybin = np.digitize(y,bins)
# 	# print('xbin',xbin)

# 	outflag = np.zeros(len(x),dtype=bool)
# 	bad_stars = []

# 	for i, xb in enumerate(bins):
# 		for j, yb in enumerate(bins):
# 			idx = np.nonzero((xbin == i) & (ybin == j))[0]

# 			stdx = np.std(dx[idx])
# 			stdy = np.std(dy[idx])
# 			stdr = np.sqrt(stdx**2 + stdy**2)

# 			meanx = np.mean(dx[idx])
# 			meany = np.mean(dy[idx])

# 			r = np.sqrt((dx[idx]-meanx)**2 + (dy[idx]-meany)**2)

# 			rx = abs(dx[idx]-meanx)
# 			ry = abs(dy[idx]-meany)
# 			out = np.nonzero((rx > 3*stdx) | (ry > 3*stdy))[0]
# 			outflag[idx[out]] = 1
# 			# gaia_id2 = gaia_id[idx] 				
# 			bad_stars.extend(gaia_id[idx[out]])

# 			# print('len rx',len(rx))
# 			# print('len out',len(out))
# 			# print('len idx',len(idx))

			
# 			# print('idx',idx)
# 			# # print('stdr', stdr, stdr.shape)
# 			# # print('r', r, r.shape)
# 			# print('where', out)
# 			# print('out',outflag[out])
	
# 	return outflag, bad_stars

def find_outliers(x,y,xref,yref,gaia_id,verbose=False):
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
					bad_stars.extend(gaia_id[include][idx[out]])

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
				bad_stars.extend(gaia_id[idx[out]])

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

	return outflag, bad_stars, m_distances	

def mahalanobis(dx,dy):
	mu = np.array([np.mean(dx),np.mean(dy)])
	x = np.column_stack((dx,dy))
	S = np.cov(x.T)
	SI = np.linalg.inv(S)
	D = np.diag(np.sqrt(np.dot(np.dot((x-mu), SI), (x-mu).T)))
	return D

def plot_quiver_all(x,y,xref,yref,include,outliers,plotdir):
	quiv_scale=200
	quiv_label = '10 pix'
	quiv_label_val = 10.0
	# plt.close('all')
	plt.figure(num=1,figsize=(6,6),clear=True)
	q = plt.quiver(x[include],y[include],(xref[include]-x[include]),(yref[include]-y[include]), color='black', scale=quiv_scale, angles='xy',width=0.0004)
	q = plt.quiver(x[outliers],y[outliers],(xref[outliers]-x[outliers]),(yref[outliers]-y[outliers]), color='red', scale=quiv_scale, angles='xy',width=0.0004)
	# quiv_label = '10 pix'
	# quiv_label_val = 10.0
	plt.quiverkey(q, 0.5, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')
	plt.xlim(-400,2448)
	plt.ylim(-400,2448)
	# plt.axis('equal')
	# plt.set_aspect('equal','box')
	plt.title('All matched stars')
	plt.savefig(plotdir + 'quiver_all.jpg', bbox_inches='tight',dpi=200)
	# plt.close()




# def load_legendre_iteration(solution_year,ref_iteration):
# 	# coefficients_2021 = '../iterating_BA/results_stacked/{}/distortion_coefficients_{}_{}.txt'.format(reference_instrument,'2021',ref_iteration)
# 	# coefficients_2020 = '../iterating_BA/results_stacked/{}/distortion_coefficients_{}_{}.txt'.format(reference_instrument,'2020',ref_iteration)
# 	# coefficients_2021 = '{}distortion_coefficients_{}_{}.txt'.format(resultDir,'2021',ref_iteration)
# 	# coefficients_2020 = '{}distortion_coefficients_{}_{}.txt'.format(resultDir,'2020',ref_iteration)



# 	# if os.path.exists(coefficients_2021):
# 	# 	print('Reading in coefficients for ref_iteration {}'.format(ref_iteration))

# 	# 	if night in ['n1','n2','n3','n4']:
# 	# 		transform = ascii.read(coefficients_2020,format='fixed_width')

# 	# 	elif night in ['n5']:
# 	# 		transform = ascii.read(coefficients_2021,format='fixed_width')
# 	# 	else:
# 	# 		print('Unknown night', night)


# 	# 	px = transform['px_val']
# 	# 	py = transform['py_val']

# 	# 	# order = 3
# 	# 	#number of coefficients in px = (order+1)^2 
# 	# 	#order = sqrt(num_coefficents) -1
# 	# 	order = int(np.sqrt(len(px))-1)
# 	coefficients = '{}distortion_coefficients_{}_{}.txt'.format(resultDir,solution_year,ref_iteration)

# 	if os.path.exists(coefficients):
# 		transform = ascii.read(coefficients,format='fixed_width')
# 		px = transform['px_val']
# 		py = transform['py_val']

# 		# order = 3
# 		#number of coefficients in px = (order+1)^2 
# 		#order = sqrt(num_coefficents) -1
# 		order = int(np.sqrt(len(px))-1)

# 	else:
# 		print('Coefficients not found for ref_iteration {}'.format(ref_iteration))
# 		return None

# 	x_domain = [0,2048] #This is not the default domain, and could mess things up?   this function doesn't get called anyway
# 	y_domain = [0,2048]
# 	tform = transforms.LegTransform(order, px, py, x_domain, y_domain, astropy_order=True)
# 	return tform



# def load_legendre(night):
# 	make solution year a parameter here.

# 	# coefficients_2021 = '../results_stacked/{}/distortion_coefficients_{}.txt'.format(reference_instrument,'2021')
# 	# coefficients_2020 = '../results_stacked/{}/distortion_coefficients_{}.txt'.format(reference_instrument,'2020')
# 	coefficients_2021 = '{}distortion_coefficients_{}.txt'.format(resultDir,'2021')
# 	coefficients_2020 = '{}distortion_coefficients_{}.txt'.format(resultDir,'2020')

# 	if os.path.exists(coefficients_2021):
# 		print('Reading in coefficients')


# 		if night in ['n1','n2','n3','n4']:
# 			transform = ascii.read(coefficients_2020,format='fixed_width')

# 		elif night in ['n5']:
# 			transform = ascii.read(coefficients_2021,format='fixed_width')
# 		else:
# 			print('Unknown night', night)


# 		px = transform['px_val']
# 		py = transform['py_val']

# 		# order = 3
# 		#number of coefficients in px = (order+1)^2 
# 		#order = sqrt(num_coefficents) -1
# 		order = int(np.sqrt(len(px))-1)
# 		# print('order',order)


# 	else:  #hardcoded if file not found.
# 		print('Using hardcoded coefficients')
# 		px = [
# 		    -15.900759376126523,
# 		     1.0364201605026848,
# 		-4.2480132859604875e-05,
# 		 2.6583671187161362e-08,
# 		  -7.96665871102779e-12,
# 		  8.356912384270864e-16,
# 		    0.02959738780768572,
# 		-0.00017377429602359357,
# 		  3.504810948108154e-07,
# 		-2.4251845608137213e-10,
# 		   6.93655490203978e-14,
# 		  -7.03160443545195e-18,
# 		-2.8355145263151783e-05,
# 		 3.4501648194713185e-07,
# 		 -6.605102285833976e-10,
# 		  4.422211586326196e-13,
# 		-1.2369493213340534e-16,
# 		 1.2339161657904676e-20,
# 		 1.6077420863164148e-08,
# 		-2.3232456244768562e-10,
# 		 4.3030662163902005e-13,
# 		  -2.80882165689627e-16,
# 		  7.704241659585834e-20,
# 		 -7.567274370242286e-24,
# 		-4.5299344365451646e-12,
# 		  6.222382681424659e-14,
# 		-1.1260007169897593e-16,
# 		  7.191259140641926e-20,
# 		 -1.935046735500897e-23,
# 		  1.868942420470575e-27,
# 		 4.1582959056782377e-16,
# 		-5.6688371394574754e-18,
# 		 1.0083414207363985e-20,
# 		 -6.304217181293066e-24,
# 		 1.6607202125669892e-27,
# 		-1.5711520442042872e-31,
# 		 ]
# 		py = [
# 		    -0.4754184802274096, 
# 		  0.0044104277258321155, 
# 		-1.1650195060774284e-05, 
# 		  8.748020034922065e-09, 
# 		 -2.689128898754918e-12, 
# 		  2.886425555007885e-16, 
# 		     0.9930218501704915, 
# 		 3.0115877014109493e-06, 
# 		   4.70148750509497e-08, 
# 		 -5.059469473711544e-11, 
# 		 1.7988616236463992e-14, 
# 		 -2.094402311176547e-18, 
# 		 2.6712941414891383e-05, 
# 		 -8.919910944339552e-08, 
# 		  5.058674840427255e-11, 
# 		  2.337261217270854e-15, 
# 		 -7.316137910387337e-18, 
# 		 1.2520194041002506e-21, 
# 		-2.4227490012256408e-08, 
# 		  9.754938238137539e-11, 
# 		 -9.326231769994538e-14, 
# 		  3.666889372703311e-17, 
# 		-5.6996379506106044e-21, 
# 		  2.073467688485445e-25, 
# 		  8.514426157741642e-12, 
# 		-3.6523749719226827e-14, 
# 		 3.9422370324078333e-17, 
# 		-1.8234552053545934e-20, 
# 		  3.669447946348159e-24, 
# 		-2.5173060076686653e-28, 
# 		-1.0253314880712769e-15, 
# 		  4.571461356122758e-18, 
# 		  -5.22524702494984e-21, 
# 		 2.5571919039766334e-24, 
# 		 -5.492974764006351e-28, 
# 		 4.1500255549075803e-32, 
# 		  ]
# 		order = 5

# 	x_domain = [0,2048]
# 	y_domain = [0,2048]
# 	tform = transforms.LegTransform(order, px, py, x_domain, y_domain, astropy_order=True)
# 	return tform

# for n in ['n1', 'n2', 'n3', 'n4']:
# 	main(n)

# main('n1')
# main('n2')
# main('n3')
# main('n4')
# main('n5')

# main('n5',fitmode='TT')   #fails to find matches.
# main('n5',fitmode='FF')
# main('n5',fitmode='FT').  #fails to find matches in ref_iteration ~3
# main('n5',fitmode='TF')

# main('2020',order=5)
# main('2021',order=5)
# for order in [4,5,6,7,8,9]:
# 	main('2020',order=order)
# 	main('2021',order=order)
for order in [2,3]:
	main('2020',order=order)
	main('2021',order=order)
