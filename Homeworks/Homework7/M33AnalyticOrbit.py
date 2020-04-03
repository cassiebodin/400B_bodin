#!/usr/bin/env python
# coding: utf-8

# # Homework 7
# ## Cassandra Bodin
# 
# 
# Template provided by: Rixin Li & G . Besla

# ## M33AnalyticOrbit Code:

# In[46]:


# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex
get_ipython().run_line_magic('matplotlib', 'inline')


# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass2 import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass


# In[60]:


class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename, ptype=2):
    # Initialize the instance of this Class with the following properties:
    # Note: ptype=2 refers to the disk particles        

        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        self.filename = filename
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        COM_M33=CenterOfMass("M33_000.txt",2)
        
        #COM position and velocity for M33
        # delta = 0.1, VolDec=4
        COMP_M33=COM_M33.COM_P(0.1,4)
        COMV_M33=COM_M33.COM_V(COMP_M33[0],COMP_M33[1],COMP_M33[2])      
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        COM_M31=CenterOfMass("M31_000.txt",2)
        
        #COM position and velocity for M33
        # delta = 0.1, VolDec=2
        COMP_M31=COM_M31.COM_P(0.1,2)
        COMV_M31=COM_M31.COM_V(COMP_M31[0],COMP_M31[1],COMP_M31[2]) 
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0= (COMP_M33 - COMP_M31).value 
        self.v0= (COMV_M33 - COMV_M31).value
        
        #store the scale lengths and masses for each component of M31:
        #for the halo (ptype=1)
        self.rhalo=63 #scale length calculated in HW5 (kpc) 
        self.Mhalo=ComponentMass("M31_000.txt", 1)*1e12 #(Msun)
        #for the disk (ptype=2)
        self.rdisk=5 #(kpc) 
        self.Mdisk=ComponentMass("M31_000.txt", 2)*1e12 #(Msun)
        #for the bulge (ptype=3)
        self.rbulge=1 #(kpc) 
        self.Mbulge=ComponentMass("M31_000.txt", 3)*1e12 #(Msun)
        
    #define function to take into account halo and bulge acceleration using Hernquist profile
    def HernquistAccel(self,M,r_a,r):  
        #Inputs:
        #M, will be either the total halo mass (Mhalo) or total bulge mass (Mbulge)
        #r_a, corresponding scale length (kpc)
        #r, a 3D relative position vector
        
        #Returns:
        #Hernquist acceleration
        
        ### **** Store the magnitude of the position vector
        rmag = np.sqrt(r[0]**2+r[1]**2+r[2]**2)
        
        ### *** Store the Acceleration
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(r_a + rmag)**2) * r --> where the last r is a VECTOR 
        Hern =  -self.G*M/(rmag *(r_a + rmag)**2) * r
        
        
        return Hern
    
    
    #define function to take into account diskacceleration using Miyamoto Nagai 1975 profile
    def MiyamotoNagaiAccel(self,M,rd,r):
        # it is easiest if you take as an input a position VECTOR  r 
        #Inputs:
        #M,the mass of the disk self.Mdisk (Msun))
        #rd, scale length of the disk self.rdisk (kpc)
        #r, a 3D relative position vector
        
        #Returns:
        #Miyamoto Nagai acceleration

        #set vriables
        zd=self.rdisk/5.0 #scale radius of disk
        R=np.sqrt(r[0]**2+r[1]**2)
        B=rd+np.sqrt(r[2]**2+zd**2)
        
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        ZSTUFF = 1/np.sqrt(r[2]**2+zd**2)
        
        a=-self.G*M/(R**2+B**2)**1.5*r*np.array([1,1,ZSTUFF])
       
        return a
        # the np.array allows for a different value for the z component of the acceleration
     
    #define a function that sums all acceleration vectors from ea galaxy component
    def M31Accel(self,r): 
        #Inputs:
        #r, a 3D relative position vector
        
        #Returns:
        #the sum of each components of acceleration based on the halo, bulge, and disk
        
        
        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        SUM=self.HernquistAccel(self.Mhalo,self.rhalo,r)+self.HernquistAccel(self.Mbulge,self.rbulge,r) +self.MiyamotoNagaiAccel(self.Mdisk,self.rdisk,r)    
        
        # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return SUM
    
    
    #define a function to integrate the acceleration vector over some amount of time
    def LeapFrog(self,dt,r,v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        #Inputs:
        #dt, timeinterval for integration
        #r, starting 3D position vector (M33-M31)
        #v, starting 3D velocity vector (M33-M31)
        
        #Returns:
        #the new position and velocity vectors
        
        # predict the position at the next half timestep
        rhalf = r+v*dt/2
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v+ self.M31Accel(rhalf)*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf+0.5*vnew*dt
        
        return rnew, vnew
    
    
    #define a function to solve equations of motion and compute the future orbit of M33 for 10Gyr in the future
    def OrbitIntegration(self, t0, dt, tmax):
        #Inputs:
        #t0, starting time
        #dt, interval of time
        #tmax, final time
        
        #Returns:
        #an array of positions and velocities written to a file

        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros((int(tmax/dt)+2,7))
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t<tmax):  # as long as t has not exceeded the maximal time 
            
            # **** advance the time by one timestep, dt
            t+=dt
           
            # **** store the new time in the first column of the ith row
            orbit[i,0] = t
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
                    
            
            rnew,vnew = self.LeapFrog(dt,orbit[i-1,1:4],orbit[i-1,4:7])
         
    
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            orbit[i,1:4]=rnew
            
            
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            orbit[i,4:7]=vnew
            
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i+=1
        
        
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function


# ## Plotting Code:

# In[61]:


#create an output file
AnOrbit=M33AnalyticOrbit("IntOrbit.txt")


# In[63]:


#perform the integration
#start:0Gyr
#end:10Gyr

#time increments: 0.5Gyr
#AnOrbit.OrbitIntegration(0, 0.5, 10)

#time increments: 0.01Gyr
AnOrbit.OrbitIntegration(0, 0.01, 10)


# In[77]:


#pull in information from IntOrbit file
intdata=np.genfromtxt('IntOrbit.txt',dtype=None,names=True)

#use data to find t, r, and v
#t is the time at each point
#sep is the seperation
#relvel is the relative velocity
t=intdata['t']
sep=np.sqrt(intdata['x']**2+intdata['y']**2+intdata['z']**2)
relvel=np.sqrt(intdata['vx']**2+intdata['vy']**2+intdata['vz']**2)


# In[78]:


#read in data files focusing on the time (columns 0)
#read in data files focusing on the positions x,y,z (columns 1,2,3)
#read in data files focusing on the velocities vx,vy,vz (columns 4,5,6)
M31_t=np.genfromtxt('Orbit_M31.txt',dtype=None,names=None,usecols=(0))
M31_pos=np.genfromtxt('Orbit_M31.txt',dtype=None,names=None,usecols=(1,2,3))
M31_vel=np.genfromtxt('Orbit_M31.txt',dtype=None,names=None,usecols=(4,5,6))

#read in data files focusing on the time (columns 0)
#read in data files focusing on the positions x,y,z (columns 1,2,3)
#read in data files focusing on the velocities vx,vy,vz (columns 4,5,6)
M33_t=np.genfromtxt('Orbit_M33.txt',dtype=None,names=None,usecols=(0))
M33_pos=np.genfromtxt('Orbit_M33.txt',dtype=None,names=None,usecols=(1,2,3))
M33_vel=np.genfromtxt('Orbit_M33.txt',dtype=None,names=None,usecols=(4,5,6))


# In[79]:


#copied function from homework 6
# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  
def MagDiff(array1, array2):
    #Inputs:
    #2 arrays of vectors 
    #    each component of the array of vectors will be subtracted
    
    #Returns:
    #the magnitude of the difference between the two vectors
    
    #Calculate the differences of each component of the two vectors/arrays
    xcomp= array1[:,0] - array2[:,0]
    ycomp= array1[:,1] - array2[:,1]
    zcomp= array1[:,2] - array2[:,2]
    
    #Calculate the magnitude of the differences of the vectors
    Magnitude=np.sqrt(xcomp**2+ycomp**2+zcomp**2)
    
    return Magnitude

# Determine the magnitude of the relative position (pos) and velocities (vel) 

# of M33 and M31
posmag_M33M31 = MagDiff(M33_pos,M31_pos)
velmag_M33M31 = MagDiff(M33_vel,M31_vel)


# In[81]:


# Plot the Orbit of the galaxies 
#################################


#initialize plot
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)



#plot the position separation between M31 and M33 as a function of time
#time is the same so can use M31_t for the time for each

plt.plot(M31_t, posmag_M33M31, color='blue', linewidth=3, label='M31-M33 simulation')
plt.plot(t, r, color='red', linewidth=3, label='M31-M33 analytic')


#set plot title
plt.title('Magnitude of Seperation Between Galaxies')

#set axes labels
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Distance (kpc)', fontsize=22)

#set axes limits
plt.xlim(0,12)
plt.ylim(0,550)

#Create and Customize Label
ax.legend(loc='upper right', fontsize='x-large')

# Save to a file
plt.savefig('MagSep.pdf', rasterized=True, dpi=350)


# In[85]:



# Plot the orbital velocities of the galaxies 
#################################

#initialize plot
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

#plot the relative velocities between M31 and M33 as a function of time
#time is the same so can use M31_t for the time for each

plt.plot(M31_t, velmag_M33M31, color='blue', linewidth=3, label='M31-M33 simulation')
plt.plot(t, relvel, color='red', linewidth=3, label='M31-M33 analytic')

#set plot title
plt.title('Magnitude of Relative Velocities Between Galaxies')

#set axes labels
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel(' Relative Velocity (km/s)', fontsize=22)

#set axes limits
plt.xlim(0,12)
plt.ylim(50,400)

#Create and Customize Label
ax.legend(loc='upper right', fontsize='x-large')

# Save to a file
plt.savefig('MagRelVel.pdf', rasterized=True, dpi=350)


# ## Questions
# 
# ### 2) How do the plots compare?
# Both the plots initially are pretty close until they begin diverging at around 1-1.25Gyr. After that point you only see one very long orbit in the analytic data, whereas there are many shorter orbits in the simulation data.
# 
# ### 3) What missing physics could make a difference?
# The analytic solution doesn't take into account the effects of the Milky Way on the M33-M31 interactions. The Milky Way is significantly larger than M33 and similar size to M31 so it would have drastic effects on the outcome of their orbits. Taking the Milky Wayinto account would also change this from a 2body orbit problem to a 3 body orbit problem.
# 
# ### 4)The MW is missing in these calculations. How might you include its effects?
# To account for the milky way we would have to alter the code to become a 3body orbital problem which takes into account each bodies interactions with each of the other bodies. New COM positions and velocites would have to be found for 3 bodies instead of 2 and then that data used to map out the seperations and relative velocities between M31 and M33. This should make it look closer to the simulation data.
# 

# In[ ]:




