#generate_pcu_grid.py

#Generates the reference grid of pcu pinholes, assuming they are perfectly machined


import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.table import Table

grid_spacing = 0.5  #in mm

grid_diameter = 44  #mm

grid_radius = grid_diameter/2
# Even number of points across grid. So odd number of gaps

#  88 pinholes across

def main():

	grid_steps = np.arange(grid_spacing/2,grid_radius,grid_spacing)
	grid_steps = np.concatenate((-grid_steps[::-1],grid_steps))
	print(len(grid_steps))
	x_array = []
	y_array = []
	id_array = []
	counter=1
	for x in grid_steps:
		for y in grid_steps:
			if ((abs(x) == 3.25) & (abs(y) == 21.75)) or ((abs(y) == 3.25) & (abs(x) == 21.75)):
				#deleting some points that aren't in the pinhole mask reference.
				continue
			if np.hypot(x,y)<=grid_radius:
				x_array.append(x)
				y_array.append(y)
				id_array.append(counter)
				counter +=1
	print(grid_steps)

	id_array = np.arange(1,len(x_array)+1)
	pinhole_size = np.ones(len(x_array))*0.024
	#save x_array and y_array to a text file, like a starlist. Import a function?

	output_table = Table([id_array,x_array,y_array,pinhole_size],names=['id','x','y','size'])
	output_table.write('pinhole_grid_positions.txt',format='ascii.fixed_width',overwrite=True)

	theta = np.linspace(0,2*math.pi,1000)
	x_circ = grid_radius*np.cos(theta)
	y_circ = grid_radius*np.sin(theta)
	
	plt.close('all')
	plt.figure(figsize=(6,6))
	plt.scatter(x_array,y_array,marker='.')
	plt.axis('equal')
	plt.show()





main()