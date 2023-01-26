import numpy as np
import astropy.units as u
from ReadFile import Read

def magnitude(qx, qy, qz):
	"""This function gives the magnitude of a vector quantity given its three components.

	Parameters
	----------
	qx, qy, qz: float
		The three components of the given quantity.

	Returns
	-------
	mag: float
		The magnitude of the vector quantity.
	"""

	mag = np.sqrt(qx**2 + qy**2 + qz**2)		# Calculating the magnitude of the quantity.

	return mag

def ParticleInfo(filename, part_type, part_no):
	"""This function extracts the following three properties of any particle from the data file given its type and number: net distance, net speed, and mass.

	Parameters
	----------
	filename: str
		Path and name of the data file.

	part_type: int
		Integer representing the particle type.
		1=Dark Matter, 2=Disk Star, 3=Bulge Star.

	part_no: int
		Number of the particle within the given particle type.

	Returns
	-------
	dist: Quantity
		Magnitude of the distance in 'kpc'.

	speed: Quantity
		Magnitude of the velocity in 'km/s'.

	mass: Quantity
		Mass in units of 'M_sun'.
	"""

	time, total_part, data = Read(filename)		# Calling Read() from ReadFile.py to read the data from the dat file.

	index = np.where(data['type'] == part_type)	# Extracting indexes of the given type of particles.
	data_new = data[index]				# Constructing a new array of particle data for given type of particles.

	x, y, z = data_new['x'][part_no-1], data_new['y'][part_no-1], data_new['z'][part_no-1]		# Finding positions of the given particle.
	x, y, z = float(x), float(y), float(z)								# Converting positions to float from str.
	vx, vy, vz = data_new['vx'][part_no-1], data_new['vy'][part_no-1], data_new['vz'][part_no-1]	# Finding velocities of the given particle.
	vx, vy, vz = float(vx), float(vy), float(vz)							# Converting velocities to float from str.

	dist = np.round(magnitude(x, y, z), 3)				# Calculating the magnitude of the distance and rounding to 3 decimal places.
	speed = np.round(magnitude(vx, vy, vz), 3)			# Calculating the magnitude of the velocity and rounding to 3 decimal places.
	mass = np.round(float(data_new['m'][part_no-1])*1e+10, 3)	# Finding mass of the given particle and rounding to 3 decimal places.

	dist = dist*u.kpc				# Attaching the kpc units.
	speed = speed*u.km/u.s				# Attaching the km/s units.
	mass = mass*u.M_sun				# Attaching the M_sun units.

	return dist, speed, mass

# Testing:
#dist, speed, mass = ParticleInfo('/home/jaymotka/400B_2023_motka/Data/MW_000.txt', 2, 100)
#print('The properties of the 100th disk partcle of the Milky Way at SnapNumber 0 is...\n')
#print('Ans. 1) 3D distance = ' + str(dist) + '.')
#print('Ans. 2) 3D velocity = ' + str(speed) + '.')
#print('Ans. 3) Mass = ' + str(mass) + '.')
#print('Ans. 4) 3D distance in lightyears = ' + str(np.round(dist.to(u.lyr), 3)) + '.\n')
