#!/usr/bin/env python
# coding: utf-8

# This program will read in data from the readfile which will return the properties for any given particle of any given type.

# In[1]:


#import libraries
#use abreviations to simplify later code
import numpy as np
import astropy.units as u
from ReadFile import Read
from astropy.constants import kpc


# In[14]:


#define a function that will extract the particle information from the ReadFile program we created.
def ParticleInfo(filename,par_type, par_num):
    #Inputs:
    #filename, we will be using MW_000.txt as our file
    #par_type, the particle type. Will use to index the information
    #par_num, the number of particles
    
    #Returns:
    #R, the calculated 3D distance using the new x,y, and z data
    #V, the calculated 3D velocity using the new vx,vy, and vz data
    #M, the calculated mass
    
    #pull data from the ReadFile using our define import
    time, tot_part, data = Read(filename)
    
    #creates an index for all the particles given a property
    index = np.where(data['type'] == par_type)
    
    #stores those components in the index as new variables 
    xnew=data['x'][index][par_num]
    ynew=data['y'][index][par_num]
    znew=data['z'][index][par_num]
    vxnew=data['vx'][index][par_num]
    vynew=data['vy'][index][par_num]
    vznew=data['vz'][index][par_num]
    mnew=data['m'][index][par_num]
    
    
    #calculate values using the new variables:
    #    R, the calculated 3D distance using the new x,y, and z data
    #    V, the calculated 3D velocity using the new vx,vy, and vz data
    #    M, the calculated mass    
    R = np.sqrt(xnew**2+ynew**2+znew**2)*u.kpc
    V = np.sqrt(vxnew**2+vynew**2+vznew**2) *u.kilometer/u.second
    M = mnew*u.solMass
    
    
    return R, V, M

#pull the values from the function and redefine them
# with defined filename=MW_000.txt, par_type=2 and par_num=99
#R is now r
#V is now v
#m is now m
r, v, m = ParticleInfo("MW_000.txt",2, 99)

#print statements (all values rounded 20 3 decimals)
#print out the radius in kpc and convert to Lyr
print("Radius is ", np.round(r,3))
print("Converted Radius is ", np.round(r.to_value(u.lyr)*u.lyr,3))

#print the velocity in km/s
print("Velocity is ",np.round(v,3))

#print mass in solar masses
#mass in the file was labeled as in 10^10 Msun so multiple by 1e10
print("Mass is ", m*1e10) 


# In[ ]:




