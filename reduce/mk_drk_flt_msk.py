
import numpy as np
import os, sys
import glob

from kai.reduce import calib
from kai.reduce import sky
# from kai.reduce import data
# from kai.reduce import util
# from kai.reduce import dar
from kai.reduce import kai_util
from kai import instruments

##########
# Change the epoch, instrument, and distortion solution.
##########
osiris = instruments.OSIRIS()

##########
# Make electronic logs
#    - run this first thing for a new observing run.
##########
def makelog_and_prep_images():
    """Make an electronic log from all the files in the ../raw/ directory.
    The file will be called nirc2.log and stored in the same directory.

    @author Jessica Lu
    @author Sylvana Yelda
    """
    kai_util.makelog('../raw', instrument=osiris)
    
    # If you are reducing OSIRIS, you need to flip the images first. 
    raw_files = glob.glob('../raw/*.fits')
    # raw_files = glob.glob('../raw/dark_frame/*.fits')
    osiris.flip_images(raw_files)

    # Download weather data we will need.
    # dar.get_atm_conditions('2023')

    return

def go_calib():
    """Do the calibration reduction.

    @author Jessica Lu
    @author Sylvana Yelda
    """

    ####################
    #
    # Calibration files:
    #     everything created under calib/
    #
    ####################
    # Darks - created in subdir darks/
    #  - darks needed to make bad pixel mask
    #  - store the resulting dark in the file name that indicates the
    #    integration time (2.8s) and the coadds (10ca).
    #    -- If you use the OSIRIS image, you must include the full filename in the list. 
    
    darkFiles = ['i230413_a003{0:03d}_flip'.format(ii) for ii in range(30, 30+1)]
    darkFiles += ['i230413_a003{0:03d}_flip'.format(ii) for ii in range(35, 38+1)]
    calib.makedark(darkFiles, 'dark_60s_1ca_1rd.fits', raw_dir=None, instrument=osiris)
# 
    # Flats - created in subdir flats
    #Kp filter
    offFiles = ['i230406_a001{0:03d}_flip'.format(ii) for ii in range(2, 2+1, 1)]
    onFiles  = ['i230406_a004{0:03d}_flip'.format(ii) for ii in range(2, 2+1, 1)]
    calib.makeflat(onFiles, offFiles, 'flat_kp_tdOpen.fits', raw_dir=None, instrument=osiris)

    # Masks (assumes files were created under calib/darks/ and calib/flats/)
    calib.makemask('dark_60s_1ca_1rd.fits', 'flat_kp_tdOpen.fits', 'supermask.fits', instrument=osiris)


makelog_and_prep_images()
go_calib()
