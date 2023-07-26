
from kai.reduce import analysis
from kai import instruments
# import pdb


class pcu_analysis(analysis.Analysis):
    def __init__(self, epoch, filt, rootDir='u/mfreeman/work/PCU/script_testing/day3/',
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
        self.starlist = self.rootDir + 'pcu.list'

        # Set up some extra starfinder keywords to optimize PSF handling.
        self.stf_extra_args = ', psfSize=2.0, trimfake=0' # 2 arcsec
        self.corrMain = 0.7
        self.corrSub = 0.5
        self.corrClean = 0.8

        ##########
        # Setup the appropriate calibration stuff.
        ##########
        self.mapFilter2Cal = {'kp': 'K', 'h': 'H', 'j': 'J', 'kn3_tdOpen': 'kn3', 'kp_tdOpen':'K', 'kn3_tdhBand': 'kn3','Hbb':'H'}
        
        # Use the default stars
        self.calStars = None

        # Choose the column based on the filter
        self.calColumn = self.mapFilter2Cal[filt]

        # Set the coo star
        self.cooStar = 'pcu'
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
    # temp.cleanFiles = ['ci211024_a019002_flip', 'ci211024_a019003_flip', 'ci211024_a019004_flip', 'ci211024_a019005_flip', 
                    # 'ci211024_a019006_flip', 'ci211024_a019007_flip', 'ci211024_a019008_flip', 'ci211024_a019009_flip', 
                    # 'ci211024_a019010_flip', 'ci211024_a019011_flip', 'ci211024_a019012_flip', 'ci211024_a019013_flip',    
                    # ] # star s1 is off the edge of 019, so using a different coo star s2
    # temp.cooStar = 's2' 
    # temp.calCooStar = temp.cooStar
    # temp.calFile = temp.rootDir + 's2_photo.dat'
    # temp.labellist = temp.rootDir + 's2_label.dat'
    # temp.starfinderCleanLoop()
    temp.starfinderClean()
    return





analyze_pcu()
