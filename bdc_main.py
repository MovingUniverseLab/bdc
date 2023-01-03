#this is the main script that calculates the distortion.

"""
Rough outline for code: follow the method in distortion_iterating.py.
But tidier.


Inputs:

Image starlists: We need to run starfinder on the images to get starlists. That could be included as a step here, or it could be a separate piece of code.
Either way, this script should take in the starlists.


Reference starlist: This will probably be the PCU grid, but could also be Hubble/Gaia catalogues. Make that an option.
Other target info? M15, M92, PCU grid

Initial transformation guesses: Very important for using the PCU. Maybe also important for Hubble/Gaia. I did end up using initial guesses. Might be hard to automate.


Several nights: The data may come in separate chunks, like with my observations being spread over several nights. Each night may have different parameters.
What should the name for these chunks be? I think nights is good, lines up well with observed nights. Maybe 'runs'? Or epochs? Epochs is broader.

Input parameters: May vary with different nights. Make a 'night' object that contains the parameters.


Create stacks: Combine images with the same PA into stacks with error measurements. Might need some intelligence here to pick the correct stacks.

Projecting velocity. Centring the solution. I think I'll stick with the 5th order polynomial, but it must be able to vary.
Maybe include functions for f_tests etc to compare.


For each night, have a text file with the parameters in it. Maybe lists of starlists too?

This may be used for Nirc2 as well, so don't hard code OSIRIS images.


I definitely need flystar.
Do I need KAI? instrument object do not include pixel dimensions.
I do need KAI for DAR and instruments.




Include some plotting functions.


run from the command line with an input script (config file), like maos. And an output folder.
Input script will contain links to nights.


"""
import argparse
import yaml
from astropy.table import Table
from astropy.io import ascii, fits



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
	hubbleData = fetch_hubble(config[night]['hubble_file'])
	if config['generate_reference_frame']:
		# median_residuals_b = []
		# mean_residuals_b = []
		# mean_residuals_b_squared = []
		# min_residuals_b = []
		# max_residuals_b = []
		# num_residuals_b = []
		# median_residuals_b_radec = []
		# mean_residuals_a = []
		refTable_current = {}	#holds most recent refTable for each night
		#this is the reference table that keeps getting updated by the reference section.
		refTable_current_filename = '{}refTable_current_{}'.format(resultDir,solution_year) 
		#make these into a dictionary of lists.
		if config['instrument'] == 'OSIRIS'
			osiris = instruments.OSIRIS()
		for ref_iteration in range(config['ref_iterations']):
			print('\n \n Ref Iteration {} \n'.format(ref_iteration))
			distortion_model = distortion_section(refTable_current)
			#section B
			new_refTable = reference_section(refTable_current,distortion_model)
			refTable_current = new_refTable
			#section A
		#final result here.
	else:
		distortion_model = distortion_section(refTable_current)
		#final result here.
		#outputs are printed to files.


def distortion_section(refTable_current, initial_distortion_correction = None):
	# --------------------------------------- Section B ---------------------------------------------------------------
	#fp_iteration starts here?
	#get the distortion solution, correct the data, pass that back into the fitter.
	#Do I want to save all fp_iterations, or just save the final fit? I think I want all iterations to be saved. So I can quickly re-run and plot them.

	current_distortion_correction = initial_distortion_correction
	tab1_initial = {'n1':[],'n2':[],'n3':[],'n4':[],'n5':[]} #is regenerated each ref_iteration
	# do I need the empty lists here in order to append to them?

	#if I am loading dist files then this can break if a set of fp_iterations was incomplete - it won't run the first iteration to generate it. Need to save as a file?

	for fp_iteration in fp_iterations:

		matched_star_table_name = '{}dist_measures_{}_{}_{}.txt'.format(resultDir,solution_year,ref_iteration,fp_iteration)
		matched_star_table = {}
		plate_scale_and_rotations = {}

		# if we want to generate a new matched_star_table:
		if config['generate_new_fp']:  
			obs_nights = [night for night in config]
			#-------------------------Night loop starts here-------------------------------
			for night in obs_nights:
				osiris_filenames = get_osiris_files(config[night]['stackDir'])
				sl2 = slice(config[night]['slice'][0],config[night]['slice'][1])
				osiris_filenames = osiris_filenames[sl2]
				if night not in refTable_current.keys():  #if we have not generated a new reference frame, load the starting one.
					if config['use_flystar_velocity']:
						refTable_H = prepare_hubble_for_flystar(hubbleData,config[night]['ra_field'],config[night]['dec_field'],config[night]['target'])
					else:
						starlist0 = load_osiris_file(config[night]['stackDir'] ,osiris_filenames[0])
						hubbleData_p = project_pos(hubbleData,starlist0,'hubble_'+config[night]['target'])
						refTable_H = prepare_hubble_for_flystar(hubbleData_p,config[night]['ra_field'],config[night]['dec_field'],config[night]['target'])
					refTable_current[night] = refTable_H.filled()
					with open(refTable_current_filename, 'wb') as temp:
						pickle.dump(refTable_current, temp)

				for i, filename in enumerate(osiris_filenames):
					if filename in config[night]['bad_files']:
						print('{} {} flagged as bad, skipping'.format(i,filename))
						continue
					starlist, refTable_d, PA = load_and_prepare_data(filename,refTable_current)
					try:
						transformation_guess = load_transformation_guess()
						starlist_corrected = starlist[:]
						# if current_distortion_correction is not None:
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
							tab1 = tab1_initial[night][i]
							j = 0
						else:
							j=i  # i=osiris index, j=gaia index


						ref_idx = np.where(tab1['ref_orig'] == True)[j]
						if len(ref_idx) <= 10:
							print(i, filename, 'Only', len(ref_idx),'matched, skipping')
							errorlist.append(filename[:-4])
							return
						with open(tform_file_1, 'wb') as temp:
							pickle.dump(tform, temp)

						print(i, filename, len(refTable_d), 'Reference stars,', len(starlist_corrected), 'OSIRIS stars,', len(ref_idx), 'matches')
						
						px = tform_inv[j].px
						py = tform_inv[j].py
						theta = math.atan2(px[2],px[1])
						scale = math.cos(theta) / px[1]
						plate_scale_and_rotations.setdefault('Filename',[]).append(filename)
						plate_scale_and_rotations.setdefault('Scale',[]).append(scale)
						plate_scale_and_rotations.setdefault('Rotation',[]).append(math.degrees(theta))
						plate_scale_and_rotations.setdefault('PA',[]).append(PA)
						# used_files.append(filename)
						# scales.append(scale)
						# rotations.append(math.degrees(theta))
						# PAs.append(PA)

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

				print('Mean scale =', np.mean(scales))
				# offset = np.array(rotations)-np.array(PAs)
				plate_scale_and_rotations['Difference'] = plate_scale_and_rotations['Rotation'] - plate_scale_and_rotations['PA']
				print('Mean rotation offset =', np.mean(plate_scale_and_rotations['Difference']))
				# transform_table = Table([used_files, scales, rotations, PAs, offset], names = ('File','Scale', 'Rotation','PA','Difference'))
				transform_table = Table(plate_scale_and_rotations)
				ascii.write(transform_table,'{}t_params_{}_{}.txt'.format(resultDir,ref_iteration,night),format='fixed_width', overwrite=True)
			
				# x_O = np.array(x_O)
				# y_O = np.array(y_O)
				# x_G = np.array(x_G)
				# y_G = np.array(y_G)
				# weights = np.array(weights)
				# idnums = np.array(idnums)
				# starlist_num = np.array(starlist_num)
				# night_col = np.array(night_col)
				# used_in_trans_1 = np.array(used_in_trans_1)

				# output_table = Table([x_O,y_O,x_G,y_G,weights,idnums,xe_O,ye_O,xe_G,ye_G,Ra_O,Dec_O,Ra_G,Dec_G, night_col,starlist_num,used_in_trans_1], names = ('x_OSIRIS','y_OSIRIS','x_REF','y_REF','Weight','REF_ID','xe_OSIRIS','ye_OSIRIS','xe_REF','ye_REF','Ra_OSIRIS','Dec_OSIRIS','Ra_REF','Dec_REF','Night','Frame','UiT'),)
				output_table = Table(matched_star_table) #elements are lists, not numpy arrays. Should be fine?
				ascii.write(output_table,matched_star_table,format='fixed_width', overwrite=True)

				# distortion_data = output_table
				# distortion_data = ascii.read(matched_star_table,format='fixed_width')
				#I generate this for each 4p loop and overwrite it.
				#always pass a first one in.
	
		else:
			print('b: Loading fit {}'.format(matched_star_table_name))
		
		# Then always load from file.
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


		# x = distortion_data['x_IM'].data
		# y = distortion_data['y_IM'].data
		# xref = distortion_data['x_REF'].data
		# yref = distortion_data['y_REF'].data
		# weights = distortion_data['Weight'].data
		# gaia_id = distortion_data['REF_ID'].data
		# xe = distortion_data['xe_IM'].data
		# ye = distortion_data['ye_IM'].data
		# xeref = distortion_data['xe_REF'].data
		# yeref = distortion_data['ye_REF'].data
		# frame = distortion_data['Frame_num'].data
		# night_c = distortion_data['Night'].data
		# ra = distortion_data['RA_IM'].data
		# dec = distortion_data['Dec_IM'].data
		# raref  = distortion_data['RA_REF'].data
		# decref  = distortion_data['Dec_REF'].data		

		# flystar rejects outliers and sets their weights to zero, which sets their errors to inf. So we put them back here.
		ind_r = np.where(distortion_data['Weight'] == 0.0)[0] #index of stars rejected with outlier_rejection, have their weights set to 0.0
		distortion_data['Weight'][ind_r] = 1/np.sqrt((distortion_data['xe_IM'][ind_r]*0.01)**2 + (distortion_data['ye_IM'][ind_r]*0.01)**2 + distortion_data['xe_REF'][ind_r]**2 + distortion_data['ye_REF'][ind_r]**2)
		print('{} stars rejected due to Flystar outlier_rejection'.format(len(ind_r)))

		#maybe don't declare all these variables, just use the table, and update the weights in place.

		outliers, bad_names, m_distances = find_outliers(distortion_data['x_IM'],distortion_data['y_IM'],distortion_data['x_REF'],distortion_data['y_REF'],distortion_data['REF_ID'])
		print('Mean of mahalanobis distances:', np.mean(m_distances))
		print('Standard deviation of mahalanobis distances:', np.std(m_distances))
		include = ~outliers 

		print('Fitting Legendre Polynomial')
		legendre_transformation = transforms.LegTransform.derive_transform(distortion_data['x_IM'][include], 
																			distortion_data['y_IM'][include], 
																			distortion_data['x_REF'][include], 
																			distortion_data['y_REF'][include], 
																			order, m=None, mref=None,init_gx=None, init_gy=None, weights=None, mag_trans=True
																			)
		#   Defines a bivariate legendre tranformation from x,y -> xref,yref using Legnedre polynomials as the basis.

		current_distortion_correction = legendre_transformation  #update current_distortion_correction with latest Legendre polynomial

		distortion_plot_print_save(distortion_data,current_distortion_correction,fp_iteration)



	return current_distortion_correction #returns the final distortion correction




def reference_section():
		# --------------------------------------- Section A ---------------------------------------------------------------

		if len(refTable_current) == 0:
			with open(refTable_current_filename, 'rb') as temp:
				refTable_current = pickle.load(temp) 

		#------------------------Night loop 1 starts here---------------------------

		obs_nights = [night for night in config]
		print(obs_nights)
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


			plotDir_n = plotDir + night + '/'

			osiris_filenames = get_osiris_files(config[night]['stackDir'])

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
						starlist = load_osiris_file(config[night]['config[night]['stackDir']'] ,filename)
						# plt.close('all')

						# ido = np.where(starlist['m'] < 15.5)
						# starlist = starlist[ido]
						fitsfile = config[night]['cleanDir'] + filename[:-12] + '.fits'

						PA = get_PA(fitsfile)
						# print('PA', PA)
						# starlist = brightest_n(starlist,170)
						starlist = mag_cut(starlist,0,minmag)
						if not filename in config[night]['dont_trim']:
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

#------------------------end----------------------





def find_outliers(x,y,xref,yref,gaia_id,verbose=False):
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
	print('Outliers found')
	# bad_names.sort()    #probably don't need this section on bad reference stars.
	# bad_count = Counter(bad_names)
	# all_count = Counter(gaia_id)
	# bad_stars = [k for (k,v) in bad_count.items() if v > 11]
	return outflag, bad_stars, m_distances	

def mahalanobis(dx,dy):
	mu = np.array([np.mean(dx),np.mean(dy)])
	x = np.column_stack((dx,dy))
	S = np.cov(x.T)
	SI = np.linalg.inv(S)
	D = np.diag(np.sqrt(np.dot(np.dot((x-mu), SI), (x-mu).T)))
	return D


def do_fp_iteration(current_distortion_correction):






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

def get_osiris_files(starfindDir):
	filelist = os.listdir(starfindDir)
	starfindFiles = fnmatch.filter(filelist,'ci*.lis')
	starfindFiles.sort()
	return starfindFiles

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


def load_and_prepare_data(filename):
	starlist = load_osiris_file(config[night]['stackDir'] ,filename)
	fitsfile = config[night]['cleanDir'] + filename[:-12] + '.fits'
	PA = get_PA(fitsfile)
	# starlist = brightest_n(starlist,170)
	starlist = mag_cut(starlist,0,config[night]['minmag'])
	if not filename in config[night]['dont_trim']:
		starlist = edge_cut(starlist,5)
	if len(starlist) == 0:
		print(i,filename, '0 stars remaining after edge cut, skipping image')
		errorlist.append(filename[:-4])
		continue
	# if recalculate_plate_scale:
	# 	xt, yt = correction.evaluate(starlist['x'],starlist['y'])
	# 	# print(starlist['y'])
	# 	starlist['x'] = xt 
	# 	starlist['y'] = yt 
	refTable_t = trim_gaia(refTable_current[night],filename,PA)    
	refTable_tf = refTable_t.filled()		#unmask, required for quiver plots

	refTable_d = dar.applyDAR(fitsfile, refTable_tf, plot=False, instrument=osiris, plotdir=plotDir_n + 'dar_a/'+ str(ref_iteration) + '/')
	return starlist, refTable_d, PA


def append_to_matched_star_table(matched_star_table,tab1,ref_idx,tform,tform_inv);
	j=0

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

	# setdefault: if the key exists, return it. Otherwise, set the key to [] and return in. Either case can then be extended.
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

	plot_image(Gx,Gy,Ox,Oy,fitsfile,plotDir_n + 'img/'+ str(ref_iteration) + '/',tab1['use_in_trans'][ref_idx])
	plot_image_dots(Gx,Gy,Ox,Oy,fitsfile,plotDir_n + 'img_d/'+ str(ref_iteration) + '/',tab1['use_in_trans'][ref_idx])
	plot_quiver(Ox,Oy,Gx,Gy,filename[:-4],plotDir_n + 'quiver/'+ str(ref_iteration) + '/',tab1['use_in_trans'][ref_idx])
	if show_plots:
		plt.show()

	return


def distortion_plot_print_save(distortion_data,legendre_transformation):

		x = distortion_data['x_IM'].data
		y = distortion_data['y_IM'].data
		xref = distortion_data['x_REF'].data
		yref = distortion_data['y_REF'].data
		weights = distortion_data['Weight'].data
		gaia_id = distortion_data['REF_ID'].data
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



def stack_frames():
	pass



def load_reference_list():
	#remember, there can be different lists in one run
	pass


def load_observed_lists():
	pass



main()









