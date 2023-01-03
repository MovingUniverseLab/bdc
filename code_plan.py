code_plan.py


# rough plan of the code layout.



Global:
Load config file.



def main():
	reference_data = fetch_hubble()  #load the reference frame data from file.

	if iterating reference:   #do the whole generating new reference frame thing.
		for ref_iteration in number_ref_iterations:
			do_distortion() # generate a distortion solution using current reference frame.
			do_reference() # generate a new reference frame using current distortion solution.
	else:
		do_distortion()


def do_distortion():  #generate a distortion solution given a reference frame and ... some other parameters.
	
	current_distortion_correction = None
	tab1 = {}
	#need reference lists and OSIRIS lists

	for fp_iteration in number_fp_iterations:
		matched_table = []
		for night in nights:
			load data?
			project_positions? 				#could make these bits into prepare_data()
			for frame in OSIRIS_Frames:
				mag cut 
				edge_cut
				apply_DAR()
				if fp_iteration > 0:
					apply_distortion_correction(current_distortion_correction)
				msc = align.MosaicToRef()
				msc.fit()
				if fp_it_num == 0
					tab1[night].append(msc.ref_table)
				tform = msc.trans_list
				tform_inv = msc.trans_list_inverse
				reference_in_pix = tform_inv(reference)
				table = reference_in_pix + tab1[night]['some columns']
				matched_table.append(table)  

		remove_outliers(matched_table)
		updated_distortion_correction = fit_legendre(matched_table)
		current_distortion_correction = updated_distortion_correction
		plot_print_save()

	return updated_distortion_correction.


def prepare_data(night)
	load data?
	project_positions? 				#could make these bits into prepare_data()
	mag cut 
	edge_cut
	apply_DAR()
	return osiris_data, hubble_data
