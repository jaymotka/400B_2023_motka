# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl
from ReadFile import Read

class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot

    def __init__(self, filename, ptype):
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            filename : `str`
                snapshot file
            ptype : `int; 1, 2, or 3`
                particle type to use for COM calculations
        '''
     
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities
        self.x = self.data['x'][self.index]         # Storing x-coordinate
        self.y = self.data['y'][self.index]         # Storing y-coordinate
        self.z = self.data['z'][self.index]         # Storing z-coordinate
        self.vx = self.data['vx'][self.index]       # Storing x-velocity
        self.vy = self.data['vy'][self.index]       # Storing y-velocity
        self.vz = self.data['vz'][self.index]       # Storing z-velocity


    def COMdefine(self,a,b,c,m):
        ''' Method to compute the COM of a generic vector quantity by direct weighted averaging.
        
        PARAMETERS
        ----------
        a : `float or np.ndarray of floats`
            first vector component
        b : `float or np.ndarray of floats`
            second vector component
        c : `float or np.ndarray of floats`
            third vector component
        m : `float or np.ndarray of floats`
            particle masses
        
        RETURNS
        -------
        a_com : `float`
            first component on the COM vector
        b_com : `float`
            second component on the COM vector
        c_com : `float`
            third component on the COM vector
        '''
        # write your own code to compute the generic COM 
        #using Eq. 1 in the homework instructions
        
        # x_com = summation x_i*m_i / summation m_i    ... Eq. 1
        
        # xcomponent Center of mass
        a_com = np.sum(a*m) / np.sum(m)       # Calculating x-component of CoM vector
        # ycomponent Center of mass
        b_com = np.sum(b*m) / np.sum(m)       # Calculating y-component of CoM vector
        # zcomponent Center of mass
        c_com = np.sum(c*m) / np.sum(m)       # Calculating z-component of CoM vector
        
        # return the 3 components separately
        return a_com, b_com, c_com
    
    
    def COM_P(self, delta, volDec):
        '''Method to compute the position of the center of mass of the galaxy 
        using the shrinking-sphere method.

        PARAMETERS
        ----------
        delta : `float, optional`
            error tolerance in kpc. Default is 0.1 kpc
        volDec: 'float'
            the amount by which the radius is decreased for each repeated iteration. 
        
        RETURNS
        ----------
        p_COM : `np.ndarray of astropy.Quantity'
            3-D position of the center of mass in kpc
        '''                                                                     

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        # write your own code below
        r_COM = np.sqrt(x_COM**2 + y_COM**2 + z_COM**2)      # Calculating the magnitude of the COM position vector.

        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        x_new = self.x - x_COM                               # Changing x-coordinates of the particle to CoM reference frame.
        y_new = self.y - y_COM                               # Changing y-coordinates of the particle to CoM reference frame.
        z_new = self.z - z_COM                               # Changing z-coordinates of the particle to CoM reference frame.
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)      # Calculating the magnitude of the position of the particle in 
                                                             # the CoM reference frame.

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        r_max = max(r_new)/volDec
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        change = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (change > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(r_new < r_max)          # Finding index for particles within the reduced radius.
            x2 = self.x[index2]                       # Storing original x-coordinates for particles within the reduced radius.
            y2 = self.y[index2]                       # Storing original y-coordinates for particles within the reduced radius.                      
            z2 = self.z[index2]                       # Storing original z-coordinates for particles within the reduced radius.
            m2 = self.m[index2]                       # Storing original mass for particles within the reduced radius.

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x2, y2, z2, m2)  # Calculating CoM coordinates for particles within the 
                                                                     # reduced radius.
            # compute the new 3D COM position
            # write your own code below
            r_COM2 = np.sqrt(x_COM2**2 + y_COM2**2 + z_COM2**2)      # Calculating the magnitude of the COM position vector
                                                                     # for particles within the reduced radius.

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            change = np.abs(r_COM - r_COM2)
            # uncomment the following line if you want to check this                                                                                               
            # print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : r_max, particle separations and COM                                        

            # reduce the volume by a factor of volDec again                                                                 
            r_max /= volDec
            # check this.                                                                                              
            #print ("maxR", r_max)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            x_new = self.x - x_COM2                         # Changing x-coordinates of the particle to new CoM reference frame.
            y_new = self.y - y_COM2                         # Changing y-coordinates of the particle to new CoM reference frame.
            z_new = self.z - z_COM2                         # Changing z-coordinates of the particle to new CoM reference frame.
            r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2) # Calculating the magnitude of the position of the particle in 
                                                            # the new CoM reference frame.

            # set the center of mass positions to the refined values                                                   
            x_COM = x_COM2
            y_COM = y_COM2
            z_COM = z_COM2
            r_COM = r_COM2

            # create an array (np.array) to store the COM position                                                                                                                                                       
            p_COM = np.array([x_COM, y_COM, z_COM])

        # set the correct units using astropy and round all values
        # and then return the COM positon vector
        # write your own code below
        p_COM = np.around(p_COM, decimals=2)*u.kpc          # Rounding the position vecotr components to 2 decimals and
                                                            # adding the units
        return p_COM
        
    def COM_V(self, x_COM, y_COM, z_COM):
        ''' Method to compute the center of mass velocity based on the center of mass
        position.

        PARAMETERS
        ----------
        x_COM : 'astropy quantity'
            The x component of the center of mass in kpc
        y_COM : 'astropy quantity'
            The y component of the center of mass in kpc
        z_COM : 'astropy quantity'
            The z component of the center of mass in kpc
            
        RETURNS
        -------
        v_COM : `np.ndarray of astropy.Quantity'
            3-D velocity of the center of mass in km/s
        '''
        
        # the max distance from the center that we will use to determine 
        #the center of mass velocity                   
        rv_max = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position (x_COM, y_COM, z_COM)
        # write your own code below
        xV = self.x*u.kpc - x_COM       # Storing the x-coordinates of the particles in the CoM reference frame. 
        yV = self.y*u.kpc - y_COM       # Storing the y-coordinates of the particles in the CoM reference frame.
        zV = self.z*u.kpc - z_COM       # Storing the z-coordinates of the particles in the CoM reference frame.
        rV = np.sqrt(xV**2 + yV**2 + zV**2) # Storing the magnitudes of the positions of the particles in the CoM reference 
                                            # frame.
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(rV < rv_max)    # Finding the index for the particles within 15 kpc of the CoM.
        
        # determine the velocity and mass of those particles within the mas radius
        # write your own code below
        vx_new = self.vx[indexV]         # Storing the x-velocity of the particles within 15 kpc of the CoM.
        vy_new = self.vy[indexV]         # Storing the y-velocity of the particles within 15 kpc of the CoM.
        vz_new = self.vz[indexV]         # Storing the z-velocity of the particles within 15 kpc of the CoM.
        m_new =  self.m[indexV]          # Storing the mass of the particles within 15 kpc of the CoM.
        
        # compute the center of mass velocity using those particles
        # write your own code below
        vx_COM, vy_COM, vz_COM = self.COMdefine(vx_new, vy_new, vz_new, m_new)  # Finding the velocities of the CoM.
        
        # create an array to store the COM velocity
        # write your own code below
        v_COM = np.array([vx_COM, vy_COM, vz_COM])     # Creating the CoM velocity vectors
        
        # return the COM vector
        # set the correct units usint astropy
        # round all values                                                                                        
        v_COM = np.around(v_COM, decimals=2)*u.km / u.s  # Rounding the components to 2 decimals and adding the units.
     
        return v_COM    