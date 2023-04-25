#pcu_guess.py

#testing the pcu transformation guess script.

import numpy as np
from astropy.io import ascii, fits
from flystar import transforms
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import math
import sys

fitsfiles = ['i230413_a003{0:03d}_flip.fits'.format(ii) for ii in range(3, 11+1)]
fitsfiles = ['i230413_a003{0:03d}_flip.fits'.format(ii) for ii in range(3, 29+1)]

# fitsfiles = ['i230405_a006002_flip.fits']

raw_dir = '/u/mfreeman/work/PCU/script_testing/day3/raw/'
log_file = '/u/mfreeman/work/PCU/script_testing/day3/PCU_20230413_060326.log'
log = ascii.read(log_file,format='basic')
print(log)
# raw_dir = '/u/mfreeman/work/PCU/script_testing/day1/raw/'
ref_file = '/u/mfreeman/code/python/bdc/pinhole_grid_positions.txt'
ref_list = ascii.read(ref_file,format='fixed_width')


def main():

	# print(ref_list)

	for fitsfile in fitsfiles:
		four_p = calculate_PCU_guess(raw_dir+fitsfile)

		transformed_ref = np.zeros((len(ref_list),2))
		transformed_ref[:,0], transformed_ref[:,1] = four_p.evaluate(ref_list['x'],ref_list['y'])
		# transformed_ref = four_p.evaluate(ref_list['x'],ref_list['y'])
		# print(transformed_ref)
		plt.figure(1,figsize=(8,8))
		img = fits.getdata(raw_dir+fitsfile)
		img[np.where(img<1)] = 1
		# plt.close('all')
		vmin = 200
		vmax = 500 #np.max(img) #*0.75
		norm = LogNorm(vmin, vmax)
		plt.imshow(img, cmap='Greys_r', norm=norm, origin = "lower", )
		plt.scatter(transformed_ref[:,0],transformed_ref[:,1],c='r',marker='.',edgecolors=None)
		plt.xlim([512,1536])
		plt.ylim([512,1536])
		plt.show()

def calculate_PCU_guess(filename,instrument='OSIRIS'):
	if instrument == 'OSIRIS':
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

main()