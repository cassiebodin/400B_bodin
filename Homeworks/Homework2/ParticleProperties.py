#!/usr/bin/env python
# coding: utf-8

# This program will read in data from the readfile which will return the properties for any given particle of any given type.

# In[1]:


#import libraries
#use abreviations to simplify later code
import numpy as np
import astropy.units as u
from ReadFile import Read


# In[2]:


#define a function that will extract the particle information 
#from the ReadFile program we created

def ParticleInfo(filename,par_type, par_num):
    #Inputs:
    #filename, we will be using MW_000.txt as our file
    #par_type, the particle type.
    #     Halo:1
    #     Disk:2
    #     Bulge:3
    #par_num, the particle number
    
    #Returns:
    #R, the calculated 3D distance using the new x, y, and z data
    #V, the calculated 3D velocity using the new vx, vy, and vz data
    #M, the calculated mass
    
    #pull data from the ReadFile using our define import
    time, tot_part, data = Read(filename)
    
    #creates an index for all the particles given a property
    index = np.where(data['type'] == par_type)
    
    #stores those components in the index as new variables 
    xnew=data['x'][index]*u.kpc
    ynew=data['y'][index]*u.kpc
    znew=data['z'][index]*u.kpc
    vxnew=data['vx'][index]*u.km/u.s
    vynew=data['vy'][index]*u.km/u.s
    vznew=data['vz'][index]*u.km/u.s
    #mass in the file was labeled as in 10^10 Msun so multiply by 1e10
    mnew=data['m'][index]*1e10*u.Msun 
    
    
    #calculate values using the new variables:
    #    R, the calculated 3D distance 
    #        using the new x, y, and z data
    #    V, the calculated 3D velocity 
    #        using the new vx, vy, and vz data
    #    M, the calculated mass 
    #        usinf the new m data
    R = np.sqrt(xnew[par_num-1]**2+ynew[par_num-1]**2+znew[par_num-1]**2)
    V = np.sqrt(vxnew[par_num-1]**2+vynew[par_num-1]**2+vznew[par_num-1]**2)
    M = mnew[par_num-1]
    
    #return the 3D distance, 3D velocity, and mass
    return R, V, M


# In[3]:


#pull the values from the function and redefine them
# with defined filename=MW_000.txt, par_type=2 and par_num=99
#R is now r
#V is now v
#m is now m
r, v, m = ParticleInfo("MW_000.txt",2, 100)


# In[4]:


#print statements (all values rounded to 3 decimals)
#print out the radius in kpc and convert to Lyr
print("Radius is ", np.round(r,3))
print("Converted Radius is ", np.round(r.to_value(u.lyr)*u.lyr,3))

#print the velocity in km/s (all values rounded to 3 decimals)
print("Velocity is ",np.round(v,3))

#print mass in solar masses (all values rounded to 3 decimals)
print("Mass is ", m) 

