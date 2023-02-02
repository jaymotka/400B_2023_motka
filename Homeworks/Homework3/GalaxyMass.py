from ReadFile import Read

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

	time, total_part, data = Read(filename)

	print(data[0:2])

ComponentMass('~/mystuff/ASTR_400B/400B_2023_motka/Data/MW_000.txt', 1)
