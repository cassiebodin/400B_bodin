#!/usr/bin/env python
# coding: utf-8

# # Homework 4
# ## Center of Mass Position and Velocity
# ### Cassandra Bodin
# 
# #### Note: This code was generated using a template given to us by Gurtina Besla

# ## 1-5 CenterOfMass Code

# In[82]:


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

#import the read file 
from ReadFile import Read


# In[83]:


class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot
    
    
    def __init__(self, filename, ptype):
    # Initialize the instance of this Class with the following properties:
    
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def COMdefine(self,a,b,c,m):
    # Function to compute the center of mass position or velocity generically
    
    # input: array (a,b,c) of positions or velocities and the mass
    # returns: 3 floats  (the center of mass coordinates)

        # write your own code to compute the generic COM using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        Acom = np.dot(a,m)/np.sum(m)
        # ycomponent Center of mass
        Bcom = np.dot(b,m)/np.sum(m)
        # zcomponent Center of mass
        Ccom = np.dot(c,m)/np.sum(m)
        
        return Acom, Bcom, Ccom
    
    
    def COM_P(self, delta):
    # Function to specifically return the center of mass position and velocity                                         
    
    # input:                                                                                                           
    #        particle type (1,2,3)  
    #            Halo:1
    #            Disk:2
    #            Bulge:3
    #        delta (tolerance)       
    
    # returns: One vector, with rows indicating:                                                                                                                                                                            
    #       3D coordinates of the center of mass position (kpc)                                                             

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        
        # compute the magnitude of the COM position vector.
        RCOM = np.sqrt(XCOM**2+YCOM**2+ZCOM**2)


        # iterative process to determine the center of mass
        ########################### 

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        xNew = self.x-XCOM
        yNew = self.y-YCOM
        zNew = self.z-ZCOM
        RNEW = np.sqrt(xNew**2+yNew**2+zNew**2)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        RMAX = max(RNEW)/2.0
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
            # select all particles within the reduced radius 
            #    starting from original x,y,z, m
            index2 = np.where(RNEW < RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]
            

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2, y2, z2, m2)
            
            # compute the new 3D COM position
            RCOM2 = np.sqrt(XCOM2**2+YCOM2**2+ZCOM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # uncomment the following line if you wnat to check this                                                                                               
            #print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            RMAX = RMAX/2.0
            # check this by uncommenting print statement                                                                                            
            #print ("Rmax is ", RMAX*u.kpc)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            xNew = x2-XCOM2
            yNew = y2-YCOM2
            zNew = z2-ZCOM2
            RNEW = np.sqrt(xNew**2+yNew**2+zNew**2)

            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position                                                                                                                                                       
            COMP = [XCOM, YCOM, ZCOM]

        # set the correct units usint astropy and round all values
        # and then return the COM positon vector
        
        #COMPr is the rounded COMP vector
        
        COMPr=np.round(COMP,2)*u.kpc
        
        return COMPr
    

    def COM_V(self, COMX,COMY,COMZ):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position
        xV = (self.x*u.kpc) - (COMX)
        yV = (self.y*u.kpc) - (COMY)
        zV = (self.z*u.kpc) - (COMZ)
        RV = np.sqrt(xV**2+yV**2+zV**2)
        
        # determine the index for those particles within the max radius
        indexV = np.where(RV <= RVMAX)

        # determine the velocity and mass of those particles within the mas radius
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew = self.m[indexV]
        
        # compute the center of mass velocity using those particles
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew, vynew, vznew, mnew)

        # create a vector to store the COM velocity
        # set the correct units usint astropy
        # round all values
    
        COMV = [VXCOM, VYCOM, VZCOM]
        
        #COMVr is the rounded COMV vector with units
        COMVr= np.round(COMV,2)*u.km/u.s

        # return the COM vector                                                                                        
        return COMVr


# In[84]:


# Create a Center of mass object for the MW, M31 and M33
# below is an example of using the class for MW
MWCOM = CenterOfMass("MW_000.txt", 2)
M31COM = CenterOfMass("M31_000.txt", 2)
M33COM = CenterOfMass("M33_000.txt", 2)

#par_type, the particle type
    #     Halo:1
    #     Disk:2
    #     Bulge:3
#thus the 2 is the information about the disk (which we will use in the homework)


# In[85]:


# below gives you an example of calling the class's functions
# MW:   store the position and velocity COM
#MW_COMP = MWCOM.COM_P(0.1)
#MW_COMV = MWCOM.COM_V(MW_COMP[0], MW_COMP[1], MW_COMP[2])


# ## 6 Testing Your Code

# ### 1 
# What is the COM position and velocity vector for MW, M31, M33 at Snapshot 0 using Disk Particles only (use 1kpc for tolarance)?
# 
# In practice, disk particles work the best for the COM determination. Recall that the MW COM should be close to the origin of the coordinate system (0,0,0)

# In[96]:


#Milky Way (MW) position and velocity vectors for disk particles
#at 1kpc
MW_COMP = np.round(MWCOM.COM_P(0.1),2)
MW_COMV = np.round(MWCOM.COM_V(MW_COMP[0], MW_COMP[1], MW_COMP[2]),2)

#M31 position and velocity vectors for disk particles
#at 1kpc
M31_COMP = np.round(M31COM.COM_P(0.1),2)
M31_COMV = np.round(M31COM.COM_V(M31_COMP[0], M31_COMP[1], M31_COMP[2]),2)

#M33 position and velocity vectors for disk particles
#at 1kpc
M33_COMP = np.round(M33COM.COM_P(0.1),2)
M33_COMV = np.round(M33COM.COM_V(M33_COMP[0], M33_COMP[1], M33_COMP[2]),2)


# In[97]:


#print results:

#MW
print("MW COM Position is", MW_COMP)
print("MW COM Velocity is", MW_COMV)

#M31
print("M31 COM Position is", M31_COMP)
print("M31 COM Velocity is", M31_COMV)

#M33
print("M33 COM Position is", M33_COMP)
print("M33 COM Velocity is", M33_COMV)


# #### Answer:
# 
# MW COM Position is [-0.87  2.39 -1.42] kpc
# 
# MW COM Velocity is [-0.47  3.41 -1.33] km / s
# 
# 
# M31 COM Position is [-377.66  611.43 -284.64] kpc
# 
# M31 COM Velocity is [ 72.85 -72.14  49.  ] km / s
# 
# 
# M33 COM Position is [-476.22  491.44 -412.4 ] kpc
# 
# M33 COM Velocity is [ 44.42 101.78 142.23] km / s

# ### 2
# What is the madnitude of the current separation and velocity between MW and M31?

# In[99]:


#Magnitude of position separation between MW and M31
sep_p=abs(MW_COMP-M31_COMP)
print("separation between MW and M31", np.round(np.linalg.norm(sep_p),2))

#Magnitude of velocity separation between MW and M31
sep_v=abs(MW_COMV-M31_COMV)
print("separation between MW and M31", np.round(np.linalg.norm(sep_v),2))


# #### Answer:
# separation between MW and M31 770.14
# 
# separation between MW and M31 116.69

# ### 3
# What is the magnitude of the current separation and velocity between M33 and M31?

# In[93]:


#Magnitude of position separation between M33 and M31
sep_p=abs(M33_COMP-M31_COMP)
print("separation between M33 and M31", np.round(np.linalg.norm(sep_p),2))

#Magnitude of velocity separation between M33 and M31
sep_v=abs(M33_COMV-M31_COMV)
print("separation between M33 and M31", np.round(np.linalg.norm(sep_v),2))


# #### Answer:
# separation between MW and M31 201.08
# 
# separation between MW and M31 199.37

# ### 4
# Given that M31 and MW are about to merge, why is the interative process to determine the COM important?

# #### Answer:
# The iterative process in determining the COM is important because as the galaxies move along their path the positions of all of the things in each galaxy also move due to influences from other objects, ie the other galaxy it is colliding with. Thus the COM of individual objects will constantly be changing, altering the overall COM of the galaxies. We have to iterate in order to take into account these alterations and find the COM.

# In[ ]:




