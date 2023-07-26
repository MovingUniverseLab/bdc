#Selects a group of psf stars from a clean or combo image, to use for psf fitting.

import numpy as np
from astropy import table
from datetime import date
from astropy.io import ascii
import os
from kai.reduce import prep_analysis







psfStars = []

psfStars = [
   [1812.78, 1814.51, 1],
   [1117.74, 1190.99, 1],
   [1256.99, 1814.37, 1],
   [1048.16, 913.93, 1],
   [1742.85, 914.45, 1],
   [839.94, 1398.48, 1],
   [1742.78, 706.66, 1],
   [1187.2, 1191.06, 0],
   [1117.72, 1121.83, 0],
]

if len(psfStars) == 0:
	prep_analysis.generate_list('/u/mfreeman/work/PCU/script_testing/day3/clean/pcu_Hbb/ci230413_a003003_flip.fits')
	#generate_list() options:
	# numPSFStars = 2 #Number of stars to get
	# scale = 0.020  # scale of telescope in arcsec/pixel, defaults to 0.009942
	# fwhm = 5 # full width half max to try, in pixels 
	# target = 'ob04036'  #says name of target in plot
else:
	prep_analysis.prepStarfinder('/u/mfreeman/work/PCU/script_testing/day3/','pcu',psfStars[0],psfStars,'Hbb')