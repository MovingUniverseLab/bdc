[DEFAULT]
nights 		= ['n1.txt','n2.txt']
reference 	= 'PCU'

centred					= True
use_flystar_velocity 	= False
manually_average_star_positions = True

[n1]
night_name 		= 'n1'
target 			= 'M15'
nightDir 		= '/u/mfreeman/work/d/n1/',
cleanDir 		= nightDir + 'clean/m15_kn3_tdOpen/'
starfindDir 	= cleanDir + 'starfinder/'
stackDir 		= cleanDir + 'stacks/'
bad_files 		= ['ci200804_a022007_flip_0.8_stf.lis','ci200804_a026012_flip_0.8_stf.lis','ci200804_a027003_flip_0.8_stf.lis',]
mag_limits 		= [6,16]
mag_tolerance 	= [2, 2, 2]