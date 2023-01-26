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
	time: float
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

time, total_part, data = Read('/f/Studies/ASTR 400B/400B_2023_motka/Data/MW_000.txt')

print('The time of this snapeshot is ' + str(time) + '.')
print(type(time))
print('The total no. of particles in the galaxy are ' + str(total_part) + '.')
print(type(total_part))
print('')
print('Testing data array...')
print(type(data))
print(data.names)
