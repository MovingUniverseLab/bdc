##################################################
#
# General Notes:
# -- python uses spaces to figure out the beginnings
#    and ends of functions/loops/etc. So make sure
#    to preserve spacings properly (indent). This
#    is easy to do if you use emacs with python mode
#    and color coding.
# -- You will probably need to edit almost every
#    single line of the go() function.
# -- If you need help on the individual function calls,
#    then in the pyraf prompt, import the module and
#    then print the documentation for that function:
#    --> print nirc2.nirc2log.__doc__
#    --> print range.__doc__
#
##################################################

# Import python and iraf modules
from pyraf import iraf as ir
import numpy as np
import os, sys
import glob

# Import our own custom modules
from kai.reduce import calib
from kai.reduce import sky
from kai.reduce import data
from kai.reduce import util
from kai.reduce import dar
from kai.reduce import kai_util
from kai import instruments
from kai.reduce import analysis

##########
# Change the epoch, instrument, and distortion solution.
##########
# epoch = 'd'
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
    # offFiles = ['i230406_a001{0:03d}_flip'.format(ii) for ii in range(2, 2+1, 1)]
    onFiles  = ['i230406_a004{0:03d}_flip'.format(ii) for ii in range(2, 2+1, 1)]
    offFiles = range(0,0)
    calib.makeflat(onFiles, offFiles, 'flat_Hbb.fits', raw_dir=None, instrument=osiris) #Actually Kp, but need the same for testing.

    # Masks (assumes files were created under calib/darks/ and calib/flats/)
    calib.makemask('dark_60s_1ca_1rd.fits', 'flat_Hbb.fits', 'supermask.fits', instrument=osiris)

def go():
    """
    Do the full data reduction.
    """
    ##########
    #
    # OB170095
    #
    ##########

    ##########
    # Kp-band reduction
    ##########

    # util.mkdir('kn3')
    # os.chdir('kn3')

    name = 'pcu' #Used for directory name

    sci_files = ['i230413_a003{0:03d}_flip'.format(ii) for ii in range(3, 29+1)]   
    # print(measurement_files)
    bad_files = ['i230406_a#####.fits']
    sci_files = [i for i in sci_files if i not in bad_files]
    # sci_files=['i200804_a014002_flip', 'i200804_a028004_flip']

    # sky_files = ['i200814_a014{0:03d}_flip'.format(ii) for ii in range(2, 11+1)]   #Kn3 band, from today
    # sky_files = ['212024_a014{0:03d}_flip'.format(ii) for ii in range(7, 11+1)]   #Kp band, from yesterday
    # sky_files = ['i211024_a005{0:03d}_flip'.format(ii) for ii in range(2, 13+1)]    #Using first set of science files to make sky.
    sky_files = sci_files[:]
    refSrc = [1531, 1277]   #hard coded for PCU. The location of the centre of the PCU when it is on axis at [90,185]

    # sky.makesky(sky_files, name, 'kp_tdhBand', instrument=osiris)
    sky.makesky_fromsci(sky_files, name, 'Hbb', instrument=osiris)
    data.clean(sci_files, name, 'Hbb', refSrc, refSrc, field=name, instrument=osiris, cent_box=50,ref_offset_method='pcu',check_ref_loc=False,fixDAR=False)
    # data.calcStrehl(sci_files, 'Hbb', field=name, instrument=osiris)
    #Currently refSrc is the centre of the pinohole mask, which may be off the frame. So calcStrehl doesn't work because we are passing the same coordinates.
    #Using a different point raises another problem - the coo shift is a translation, which only works when strSrc is the same as refSrc, or there is no rotation. And we expect rotation here.

    # data.combine(sci_files, 'kn3_tdhBand', epoch, field=name, trim=0, weight='strehl', submaps=3, instrument=osiris)
    # os.chdir('../')

class pcu_analysis(analysis.Analysis):
    def __init__(self, epoch, filt, rootDir='/u/mfreeman/work/d/',
                 epochDirSuffix=None, imgSuffix=None, cleanList='c.lis', alignMagCut=' -m 20 ',
                 instrument=instruments.default_inst):
        """
        For OB170095 reduction:

        epoch -- '11may' for example
        filt -- 'kp', 'lp', or 'h'
        """
        # Initialize the Analysis object
        analysis.Analysis.__init__(self, epoch, filt=filt,
                                     rootDir=rootDir,
                                     epochDirSuffix=epochDirSuffix, imgSuffix=imgSuffix,
                                     cleanList=cleanList, instrument=instrument, airopa_mode='single')

        # Use the field to set the psf starlist
        # self.starlist = self.rootDir + 'n1/s1_psf.list'
        self.starlist = self.rootDir + 'n5/Hub_psf_n5.list'

        # Set up some extra starfinder keywords to optimize PSF handling.
        self.stf_extra_args = ', psfSize=2.0, trimfake=0' # 2 arcsec
        self.corrMain = 0.7
        self.corrSub = 0.5
        self.corrClean = 0.8

        ##########
        # Setup the appropriate calibration stuff.
        ##########
        self.mapFilter2Cal = {'kp': 'K', 'h': 'H', 'j': 'J', 'kn3_tdOpen': 'kn3', 'kp_tdOpen':'K', 'kn3_tdhBand': 'kn3'}
        
        # Use the default stars
        self.calStars = None

        # Choose the column based on the filter
        self.calColumn = self.mapFilter2Cal[filt]

        # Set the coo star
        self.cooStar = 's1'
        self.calCooStar = self.cooStar

        # Override some of the default parameters
        self.calFlags = '-f 1 -R -s 1 --searchMag=2.0 '
        self.calFile = self.rootDir + 's1_photo.dat'

        self.labellist = self.rootDir + 's1_label.dat'
        self.orbitlist = None

        # Fix align flags. Otherwise, align is using too many faint stars.
        self.alignFlags = '-R 3 -v -p -a 2 ' + alignMagCut

        self.plotPosMagCut = 20.0

        self.stfFlags = 'makePsf=1, makeRes=1, makeStars=1, fixPsf=1, trimfake=1, '

        return

def analyze_pcu():
    epoch = '/u/mfreeman/work/PCU/script_testing/day3/'  #the directory containing /clean, /raw, /reduce etc.
    print('osiris.name = ', osiris.name)
    temp = pcu_analysis(epoch, 'Hbb', instrument=osiris)
    print('temp.instrument.name = ', temp.instrument.name)
    temp.cleanFiles = ['ci211024_a019002_flip', 'ci211024_a019003_flip', 'ci211024_a019004_flip', 'ci211024_a019005_flip', 
                    'ci211024_a019006_flip', 'ci211024_a019007_flip', 'ci211024_a019008_flip', 'ci211024_a019009_flip', 
                    'ci211024_a019010_flip', 'ci211024_a019011_flip', 'ci211024_a019012_flip', 'ci211024_a019013_flip',    
                    ] # star s1 is off the edge of 019, so using a different coo star s2
    temp.cooStar = 's2' 
    temp.calCooStar = temp.cooStar
    temp.calFile = temp.rootDir + 's2_photo.dat'
    temp.labellist = temp.rootDir + 's2_label.dat'
    temp.starfinderCleanLoop()
    # temp.starfinderClean()
    return

makelog_and_prep_images()
go_calib()
go()

analyze_pcu()

