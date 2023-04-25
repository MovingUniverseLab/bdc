
#pcu_focus_test.py
import os
import fnmatch
# import imaka
# from imaka.reduce import reduce_fli
import numpy as np
from astropy.io import ascii, fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize



raw_folder = '/u/mfreeman/work/PCU/script_testing/day3/raw/'

# mask_file = raw_folder + 'test_mask.fits'

mask_file = './calib/masks/' + 'supermask.fits'

# def make_mask():
# 	test_dark= raw_folder + 'i230406_a001002.fits'  #not actually a dark. AO hatch closed, but no dark filter.
# 	dark = fits.getdata(test_dark)
# 	# threshold = np.median(dark)+ 5*np.std(dark)
# 	threshold = 650
# 	# hot_pixels = np.where(dark>threshold)
# 	hot_pixels = np.where((dark>650) | (dark<100))
# 	mask_array = np.zeros([2048,2048])
# 	mask_array[hot_pixels] = 1
# 	print(np.sum(mask_array), 'hot pixels masked outside {} - {}'.format(300,650))
# 	mask_fits = fits.PrimaryHDU(data=mask_array)
# 	mask_fits.writeto(mask_file,overwrite=True)


def reduce_PCU_data():

	# osiris_filenames = get_osiris_files(raw_folder)

	# using_numbers = np.arange(4007,4056+1,1)

	# measurement_files = [filename for filename in osiris_filenames if int(filename[-11:-5]) in using_numbers]
	# measurement_files = ['i230413_a003{0:03d}_flip.fits'.format(ii) for ii in range(3, 29+1)]	
	measurement_files = ['i230413_a003{0:03d}_flip.fits'.format(ii) for ii in range(5, 5+1)]	

	# print(measurement_files)
	bad_files = ['i230406_a#####.fits']
	measurement_files = [i for i in measurement_files if i not in bad_files]
	# print(measurement_files)

	# exposure_times = 10
	# z_values = dict(zip(measurement_files,z_position))
	# exposures = dict(zip(measurement_files,exposure_times))
	fwhm_mean = []
	fwhm_std = []

	flat = fits.getdata('./calib/flats/' + 'flat_kp_tdOpen.fits')
	dark = fits.getdata('./calib/darks/' + 'dark_60s_1ca_1rd.fits')

	for filename in measurement_files:	
		# reduce_fli.find_stars_single(raw_folder+filename,fwhm=5, threshold=6, N_passes=4, plot_psf_compare=False, mask_file=mask_file, sharp_lim=0.6, peak_max=20000,round_lim=0.7,dark=dark,flat=flat)
		# find_stars_single(raw_folder+filename,fwhm=5, threshold=5, N_passes=4, plot_psf_compare=False, mask_file=mask_file, sharp_lim=0.6, peak_max=20000,round_lim=0.7,dark=dark,flat=flat)
		pass
	print('Done finding stars')

	os.makedirs('./plots/', exist_ok=True)		
	for filename in measurement_files:
		print('Saving {}'.format(filename))
		# starlist_file = test_files[0][:-5] + '_stars.txt'
		starlist = ascii.read(raw_folder+filename[:-5]+'_stars.txt',format='basic')
		img = fits.getdata(raw_folder+filename)
		hdr = fits.getheader(raw_folder+filename)
		plt.close('all')
		plt.figure(0,figsize=(8,8))
		# norm = LogNorm(1000, 10000)
		norm=Normalize(0000, 2000)
		# v_range = (0,0.2)
		# w_lim = [0,2.4e6]
		plt.imshow(img, cmap='viridis', origin='lower',norm=norm)
		plt.savefig('./plots/script_test_{}_{:.2f}_{:.2f}.jpg'.format(filename[12:15],float(hdr['PCSFX']),float(hdr['PCSFY'])),bbox_inches='tight',dpi=200)
		plt.plot(starlist['xcentroid'],starlist['ycentroid'],'.r')
		plt.savefig('./plots/script_test_{}_{:.2f}_{:.2f}marked.jpg'.format(filename[12:15],float(hdr['PCSFX']),float(hdr['PCSFY'])),bbox_inches='tight',dpi=200)
		plt.show()

	for filename in measurement_files:
		starlist = ascii.read(raw_folder+filename[:-5]+'_stars.txt',format='basic')
		# img = fits.getdata(raw_folder+filename)
		# hdr = fits.getheader(raw_folder+filename)
		# plt.close('all')
		plt.figure(1,figsize=(8,8))
		plt.scatter(starlist['xcentroid'],starlist['ycentroid'],s=1,c='r',marker='.')
	plt.savefig('./plots/all_stars.jpg'.format(),bbox_inches='tight',dpi=200)
	plt.show()


	for filename in measurement_files:
		starlist = ascii.read(raw_folder+filename[:-5]+'_stars.txt',format='basic')
		fwhms = np.mean((starlist['x_fwhm'], starlist['y_fwhm']),axis=0)
		fwhm_mean.append(np.mean(fwhms))
		fwhm_std.append(np.std(fwhms))

	# z_position = np.asarray(z_position)
	fwhm_mean = np.asarray(fwhm_mean)
	fwhm_std = np.asarray(fwhm_std)

	print('Mean FWHM = {}'.format(np.mean(fwhm_mean)))
	return


def list_positions(day):
	if day == 1:
		raw_folder = '/u/mfreeman/work/PCU/script_testing/day1/raw/'
	elif day == 2:
		raw_folder = '/u/mfreeman/work/PCU/script_testing/day2/raw/'

	osiris_filenames = get_osiris_files(raw_folder)

	print('{}\t{}\t{}'.format('Filename','PCSFX','PCSFY'))
	for filename in osiris_filenames:
		hdr = fits.getheader(raw_folder+filename)
		print('{}\t{}\t{}'.format(filename,hdr['PCSFX'],hdr['PCSFY']))

def get_osiris_files(directory):
	filelist = os.listdir(directory)
	starfindFiles = fnmatch.filter(filelist,'i*.fits')
	starfindFiles.sort()
	return starfindFiles


def make_darks(selected_files):
	darks = mean_of_images(selected_files[0:10])
	return darks


def mean_of_images(file_list):
	img = []
	for file in file_list:
		img.append(fits.getdata(raw_folder+file))
	img=np.asarray(img)
	median_dark = np.median(img,axis=0)
	return(median_dark)



#This function from IMAKA/reduce_fli.py

from photutils import DAOStarFinder
import multiprocessing as mp
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.modeling import models, fitting
import scipy.ndimage
import warnings
from astropy.utils.exceptions import AstropyWarning

def find_stars_single(img_file, fwhm, threshold, N_passes, plot_psf_compare, mask_file, sharp_lim, peak_max,round_lim,dark=None,flat=None):
    pid = mp.current_process().pid
    
    print(f'  p{pid} - Working on image: {img_file}')
    img, hdr = fits.getdata(img_file, header=True, ignore_missing_end=True)

    if dark is not None:
        print(f'  p{pid} - Subtracting dark')
        # img = subtract_dark(img_file,dark)
        img=img-dark

    if flat is not None:
        print(f'  p{pid} - Dividing by flat')
        normalized_flat = flat/np.mean(flat)
        img=img/normalized_flat
  
    mask = fits.getdata(mask_file).astype('bool')
    img = np.ma.masked_array(img, mask=mask)
    fwhm_curr = fwhm

    # Calculate the bacgkround and noise (iteratively)
    print(f'  p{pid} - Calculating background')
    bkg_threshold_above = 1
    bkg_threshold_below = 3

    good_pix = np.where(np.isfinite(img))
    
    for nn in range(5):
        bkg_mean = img[good_pix].mean()
        # bkg_mean = np.ma.median(img[good_pix])
        bkg_std = img[good_pix].std()

        bad_hi = bkg_mean + (bkg_threshold_above * bkg_std)
        bad_lo = bkg_mean - (bkg_threshold_below * bkg_std)

        good_pix = np.where((img > bad_lo) & (img < bad_hi))
    
    bkg_mean = img[good_pix].mean()
    print(np.ma.median(img))
    bkg_std = img[good_pix].std()
    img_threshold = threshold * bkg_std 
    print(f'    p{pid} - Bkg = {bkg_mean:.2f} +/- {bkg_std:.2f}')
    print(f'    p{pid} - Bkg Threshold = {img_threshold:.2f}')

    # Detect stars
    print(f'     p{pid} - Detecting Stars')

    # Each pass will have an updated fwhm for the PSF.
    for nn in range(N_passes):
        print(f'    p{pid} - Pass {nn:d} assuming FWHM = {fwhm_curr:.1f}')
        daofind = DAOStarFinder(fwhm=fwhm_curr, threshold = img_threshold, exclude_border=True, sharplo=0, sharphi=sharp_lim, peakmax=peak_max,roundlo=-round_lim,roundhi=round_lim)
        sources = daofind(img - bkg_mean, mask=mask)
        print(f'    p{pid} - {len(sources)} sources found, now fitting for FWHM.')
        
        # Calculate FWHM for each detected star.
        x_fwhm = np.zeros(len(sources), dtype=float)
        y_fwhm = np.zeros(len(sources), dtype=float)
        theta = np.zeros(len(sources), dtype=float)
        # Calculate measure of fits
        fvu = np.zeros(len(sources), dtype=float)
        lss = np.zeros(len(sources), dtype=float)
        mfr = np.zeros(len(sources), dtype=float)
    
        # We will actually be resampling the images for the Gaussian fits.
        resamp = 1 #2 #BUG - changed this value for testing bin1 open loop
        
        cutout_half_size = int(round(fwhm_curr * 3.5))
        cutout_size = 2 * cutout_half_size + 1

        # Define variables to hold final averages PSFs
        final_psf_obs = np.zeros((cutout_size*resamp, cutout_size*resamp), dtype=float)
        final_psf_mod = np.zeros((cutout_size*resamp, cutout_size*resamp), dtype=float)
        final_psf_count = 0

        # Setup our gaussian fitter with some good initial guesses.
        sigma_init_guess = fwhm_curr * gaussian_fwhm_to_sigma
        g2d_model = models.Gaussian2D(1.0, cutout_half_size*resamp, cutout_half_size*resamp,
                                          sigma_init_guess*resamp, sigma_init_guess*resamp, theta=0,
                                          bounds={'x_stddev':[0, fwhm*2*resamp], 
                                                  'y_stddev':[0, fwhm*2*resamp], 
                                                  'amplitude':[0, 2]})
        c2d_model = models.Const2D(0.0)
        
        the_model = g2d_model + c2d_model
        the_fitter = fitting.LevMarLSQFitter()
        
        cut_y, cut_x = np.mgrid[:cutout_size, :cutout_size]
        
        fails=0

        for ss in range(len(sources)):
            x_lo = int(round(sources[ss]['xcentroid'] - cutout_half_size))
            x_hi = x_lo + cutout_size
            y_lo = int(round(sources[ss]['ycentroid'] - cutout_half_size))
            y_hi = y_lo + cutout_size

            cutout_tmp = img[y_lo:y_hi, x_lo:x_hi].astype(float)
            if ((cutout_tmp.shape[0] != cutout_size) | (cutout_tmp.shape[1] != cutout_size)):
                # Edge source... fitting is no good
                continue
        
            # Oversample the image
            cutout_resamp = scipy.ndimage.zoom(cutout_tmp, resamp, order = 1)
            #cutout_resamp /= cutout_resamp.sum() #normed sum to 1
            cutout_resamp /= cutout_resamp.max() #normed peak to 1 # BUG: what if bright outlier?
            cut_y_resamp, cut_x_resamp = np.mgrid[:cutout_size*resamp, :cutout_size*resamp]

            # Fit a 2D gaussian + constant
            with warnings.catch_warnings():
                # Suppress warnings... too many.
                warnings.simplefilter("ignore", category=UserWarning)
                warnings.simplefilter("ignore", category=AstropyWarning)
                try:
                    g2d_params = the_fitter(the_model, cut_x_resamp, cut_y_resamp, cutout_resamp, 
                                            epsilon=1e-12, acc=1e-12, maxiter=300, weights=None) #added values for better fit
                except RuntimeError as err:
                    print(err)
                    print('Skipping star {}'.format(ss))
                    fails+=1
                    continue
            g2d_image = g2d_params(cut_x_resamp, cut_y_resamp)

            # Catch bad fits and ignore. 
            if (np.isnan(g2d_params.x_mean_0.value) or
                (np.abs(g2d_params.x_mean_0.value) > (cutout_size * resamp)) or
                (np.abs(g2d_params.y_mean_0.value) > (cutout_size * resamp))):
                print(f'      p{pid} - Bad fit for {ss}')
                continue

            minimum_flux_limit = 8   #Seems to be ~10 for faint PCU images.

            # Add to our average observed/model PSFs
            if sources['flux'][ss] > minimum_flux_limit:
                final_psf_count += 1
                final_psf_obs += cutout_resamp
                final_psf_mod += g2d_image
                
            # Save the FWHM and angle.
            x_fwhm[ss] = g2d_params.x_stddev_0.value / gaussian_fwhm_to_sigma / resamp
            y_fwhm[ss] = g2d_params.y_stddev_0.value / gaussian_fwhm_to_sigma / resamp
            theta[ss] = g2d_params.theta_0.value
            # calc residuals - based on FWHM
            # relevant part of cutout
            mid_ss = cutout_resamp.shape[0]/2
            x_hi_ss = int(round(mid_ss+x_fwhm[ss])); x_lo_ss = int(round(mid_ss-x_fwhm[ss]))
            y_hi_ss = int(round(mid_ss+y_fwhm[ss])); y_lo_ss = int(round(mid_ss-y_fwhm[ss]))
            cutout_resamp_cut = cutout_resamp[y_lo_ss:y_hi_ss, x_lo_ss:x_hi_ss] 
            # fit metrics
            diff_img_ss = cutout_resamp_cut - g2d_image[y_lo_ss:y_hi_ss, x_lo_ss:x_hi_ss]
            PSF_mean_ss = np.mean(cutout_resamp_cut)
            residual_ss = np.sum(diff_img_ss**2) # Least Squares Sum (LSS)
            med_fr_ss = np.median(np.abs(diff_img_ss / cutout_resamp_cut)) # median fractional residual (MFR)
            fvu_ss = residual_ss / np.sum((cutout_resamp_cut - PSF_mean_ss)**2)  # fraction of variance unexplained (FVU)
            # Save the fit
            lss[ss] = residual_ss
            fvu[ss] = fvu_ss
            mfr[ss] = med_fr_ss

            if (plot_psf_compare == True) and (x_lo > 200) and (y_lo > 200):
                #plt.figure(4, figsize=(6, 4))
                vmin = cutout_resamp.min()
                vmax = cutout_resamp.max()

                plt.figure(4, figsize=(12,3))
                plt.clf()
                # 1. Cut out Source
                plt.subplot(1,4,1)
                plt.imshow(cutout_resamp, origin='lower',
                   vmin=vmin, vmax=vmax)
                plt.gca().add_patch(Rectangle((x_lo_ss, y_lo_ss),x_hi_ss-x_lo_ss,y_hi_ss-y_lo_ss,
                    edgecolor='red',
                    facecolor='none',
                    lw=2))
                plt.colorbar(fraction=0.046, pad=0.05)
                plt.title(f'Image (resamp={resamp:d})')
                # 2. Model of source
                plt.subplot(1,4,2)
                plt.imshow(g2d_image, origin='lower',
                           vmin=vmin, vmax=vmax)
                plt.gca().add_patch(Rectangle((x_lo_ss, y_lo_ss),x_hi_ss-x_lo_ss,y_hi_ss-y_lo_ss,
                        edgecolor='red',
                        facecolor='none',
                        lw=2))
                plt.colorbar(fraction=0.046, pad=0.05)
                plt.title(f'Model (resamp={resamp:d})')
                # 3. Residual - Subtraction
                plt.subplot(1,4,3)
                plt.imshow(cutout_resamp - g2d_image, origin='lower',
                   vmin=-vmax/6, vmax=vmax/6)
                plt.gca().add_patch(Rectangle((x_lo, y_lo),x_hi-x_lo,y_hi-y_lo,
                    edgecolor='red',
                    facecolor='none',
                    lw=1))
                plt.colorbar(fraction=0.046, pad=0.04)
                plt.title(f"Data-Model (resamp={resamp:d})")
                # 4. Residual - Fraction
                plt.subplot(1,4,4)
                plt.subplots_adjust(left=0.08)
                plt.imshow((cutout_resamp - g2d_image) / cutout_resamp, vmin=-1, vmax=1) # take out outliers?
                plt.colorbar(fraction=0.046, pad=0.05)
                plt.title('Residual fraction')
                plt.suptitle(f"Source {ss} fit, FWHM x: {x_fwhm[ss]:.2f} y: {y_fwhm[ss]:.2f} | LSS {residual_ss:.2e} | FVU {fvu_ss:.2e} | MFR {med_fr_ss:.2e}")
                plt.tight_layout()
                plt.pause(0.05)
                
                pdb.set_trace()
                
            # Some occasional display
            if (plot_psf_compare == True) and (ss % 250 == 0):
                plt.figure(2, figsize=(8, 3))
                plt.clf()
                plt.subplot(1,2,1)
                plt.subplots_adjust(left=0.08)
                plt.imshow(final_psf_obs)
                plt.colorbar(fraction=0.25)
                plt.title(f'Obs PSF (resamp = {resamp:d})')
                
                plt.subplot(1,2,2)
                plt.subplots_adjust(left=0.08)
                plt.imshow(final_psf_mod)
                plt.colorbar(fraction=0.25)
                #plt.axis('equal')
                plt.title(f'Mod PSF (resamp = {resamp:d})')
                plt.suptitle(f"Observed vs. Model PSF average fit")
                plt.pause(0.05)

                print(f'    p{pid} - ss={ss} fwhm_x={x_fwhm[ss]:.1f} fwhm_y={y_fwhm[ss]:.1f}')

                
        print('Failed for {} stars'.format(fails))
        sources['x_fwhm'] = x_fwhm
        sources['y_fwhm'] = y_fwhm
        sources['theta'] = theta
        sources['LSS'] = lss
        sources['FVU'] = fvu
        sources['MFR'] = mfr

        # Save the average PSF (flux-weighted). Note we are making a slight mistake
        # here since each PSF has a different sub-pixel position... still same for both
        # obs and model
        final_psf_obs /= final_psf_count
        final_psf_mod /= final_psf_count
        final_psf_obs /= final_psf_obs.sum()
        final_psf_mod /= final_psf_mod.sum()
        # saving psf
        img_dir_name, img_file_name = os.path.split(img_file)
        psf_dir = img_dir_name + '/psf/'
        os.makedirs(psf_dir, exist_ok=True)		
        fits.writeto(psf_dir+img_file_name.replace('.fits', '_psf_obs.fits'), final_psf_obs, hdr, overwrite=True)
        fits.writeto(psf_dir+img_file_name.replace('.fits', '_psf_mod.fits'), final_psf_mod, hdr, overwrite=True)
        #TODO: make starlist specific

        # Drop sources with flux (signifiance) that isn't good enough.
        # Empirically this is <1.2
        # Also drop sources that couldn't be fit.
        good = np.where((sources['flux'] > minimum_flux_limit) & (sources['x_fwhm'] > 0) & (sources['y_fwhm'] > 0))[0]
        sources = sources[good]

        # Only use the brightest sources for calculating the mean. This is just for printing.
        idx = np.where(sources['flux'] > minimum_flux_limit)[0]
        x_fwhm_med = np.median(sources['x_fwhm'][idx])
        y_fwhm_med = np.median(sources['y_fwhm'][idx])
        
        print(f'      p{pid} - Number of sources = {len(sources)}')
        print(f'      p{pid} - Median x_fwhm = {x_fwhm_med:.1f} +/- {sources["x_fwhm"].std():.1f}')
        print(f'      p{pid} - Median y_fwhm = {y_fwhm_med:.1f} +/- {sources["y_fwhm"].std():.1f}')

        fwhm_curr = np.mean([x_fwhm_med, y_fwhm_med])

        # formats = {'xcentroid': '%8.3f', 'ycentroid': '%8.3f', 'sharpness': '%.2f',
        #            'roundness1': '%.2f', 'roundness2': '%.2f', 'peak': '%10.1f',
        #            'flux': '%10.6f', 'mag': '%6.2f', 'x_fwhm': '%5.2f', 'y_fwhm': '%5.2f',
        #            'theta': '%6.3f', 'LSS': '%5.2f', 'FVU': '%5.2f','MFR': '%5.2f',}

        sources.rename_column('id', 'name')                         #Flystar requires these column names.
        sources.rename_column('xcentroid', 'x')
        sources.rename_column('ycentroid', 'y')
        sources.rename_column('mag', 'm')
        formats = {'x': '%8.3f', 'y': '%8.3f', 'sharpness': '%.2f',
                   'roundness1': '%.2f', 'roundness2': '%.2f', 'peak': '%10.1f',
                   'flux': '%10.6f', 'm': '%6.2f', 'x_fwhm': '%5.2f', 'y_fwhm': '%5.2f',
                   'theta': '%6.3f', 'LSS': '%5.2f', 'FVU': '%5.2f','MFR': '%5.2f',}    

        sources.write(img_file.replace('.fits', '_stars.txt'), format='ascii.fixed_width',
                          delimiter=None, bookend=False, formats=formats, overwrite=True)

    

    return




if __name__ == '__main__':
	# make_mask()
	# list_positions(1)
	# list_positions(2)


	# makelog_and_prep_images()
	# go_calib()
	reduce_PCU_data()





