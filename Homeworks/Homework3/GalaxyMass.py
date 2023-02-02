from ReadFile import Read
import numpy as np
import astropy.units as u

def ComponentMass(filename, part_type):
	"""This function finds the total mass of any desired galaxy component.

        Parameters
        ----------
        filename: str
                Path and name of the data file.

        part_type: int
                Integer representing the particle type.
                1=Dark Matter, 2=Disk Star, 3=Bulge Star.

        Returns
        -------
        mass: astropy.Quantity
		Total mass of the desired galaxy component in 10^12 M_Sun.
	"""

	time, total_part, data = Read(filename)                 # Calling Read() from ReadFile.py

	index = np.where(data['type']==part_type)               # Finding indexes of the given type of particles.
	mass_data = data['m'][index]                            # Making an array of masses of the desired type of particles.

	mass = np.sum(mass_data)*1e+10*u.Msun                   # Taking the sum of masses of the desired type of particles.

	mass = np.round(mass.to(1e+12*u.Msun), decimals=3)      # Converting to the units of 10^12 M_Sun and rounding to 3 decimal digits.

	return mass


# Testing:
#table = []
#for file in ['MW_000.txt', 'M31_000.txt', 'M33_000.txt']:
#	row = [file.strip('_000.txt')]
#	total_mass = 0
#	stellar_mass = 0
#	for i in range(1,4):
#		mass = ComponentMass('/home/jaymotka/400B_2023_motka/Data/' + file, i)
#		row.append(mass)
#		total_mass += mass
#		if i != 1:
#			stellar_mass += mass
#	row.append(total_mass)
#	row.append(np.round(stellar_mass/total_mass, decimals=3))
#	table.append(row)
#table = np.array(table)
#print(table)
