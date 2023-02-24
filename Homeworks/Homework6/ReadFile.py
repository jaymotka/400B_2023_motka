import numpy as np
import astropy.units as u

def Read(filename):
	"""This function reads the given data file and returns the time of the snapshot, total no of particle in the snapshot, and particle data.

	Parameters
	----------
	filename: str
		Path and name of the data file.

	Returns
	-------
	time: Quantity
		Time of the snapshot in 'Myrs'.

	total_part: int
		Total number of particles in the galaxy snapshot.

	data: numpy.array
		The particle information as a numpy array. 
		The column labels are: 'type, m, x, y, z, vx, vy, vz'.
	"""

	file = open(filename, 'r')	# Opening the data file.

	line1 = file.readline()		# Reading the first line from the file.
	label, value = line1.split()    # Extracting the time from the first line.
	time = float(value)*u.Myr	# Converting the time from str to float and storing it in Myr units.

	line2 = file.readline()		# Reading the second line from the file.
	label, value = line2.split()	# Extracting the total no. of particles from the second line.
	total_part = int(value)		# Converting the total no. of particles from str to int.

	file.close()			# Closing the data file.

	data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)	# Extracting the particle information in a numpy array.

	return time, total_part, data

# Testing:
#time, total_part, data = Read('/home/jaymotka/400B_2023_motka/Data/MW_000.txt')
#print('Testing ReadFile.py (Question 3)... \n')
#print('The time of this snapeshot is ' + str(time) + '.')
#print('The total no. of particles in this galaxy snapshot are ' + str(total_part) + '.\n')
#print('Testing data array...')
#print('The column labels are: ' + str(data.dtype.names) + '.\n')
#print('The properties of second particles...')
#print('Type: ' + str(data['type'][1]) + ' (Dark Matter).')
#print('Mass: ' + str(data['m'][1]) + 'e+10 M_Sun.')
#print('Possitions: x=' + str(data['x'][1]) + ' kpc, y=' + str(data['y'][1]) + ' kpc, z=' + str(data['z'][1]) + ' kpc.')
#print('Velocities: vx=' + str(data['vx'][1]) + ' km/s, vy=' + str(data['vy'][1]) + ' km/s, vz=' + str(data['vz'][1]) + ' km/s.\n')
