#!/usr/bin/env python
# coding: utf-8

# This code will read a file MW_000.txt. It will print out the time, total number of particles as variables, particle type, mass, x,y,z, and vx,vy,vz.

# In[ ]:


#import libraries
#use abreviations to simplify later code
import numpy as np
import astropy.units as u


# In[ ]:


#define a function where we will read, open, and print data from the file MW_000.txt
def Read(MW_000.txt):
    
    #opens and reads file
    file = open(MW_000.txt,'r') 
    
    #reading line 1, storing the time
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*10.0*u.Myr
    
    #reading line 2, storing the total number of particles (tot_part)
    line2 = file.readline()
    label, value = line2.split()
    tot_part = float(value)*10.0*u.Myr
    
    #close the file
    file.close
    
    #store rest of file in np.genfromtxt to use the header information
    #parameters:
        #"dtype=None" splits line using white spaces
        #"names=True" creates arrays to store data with the right labels
            #labels are: type, m, x, y, z, vx, vy, vz
        #"skip_header=3" skips the first 3 lines
    data = np.genfromtxt(MW_000.txt, dtype=None, names=True, skip_header=3)
    
    return time, tot_part, data
    
print(data['type'][1])
print(data['m'][1])
print(data['x'][1])
print(data['y'][1])
print(data['z'][1])
print(data['vx'][1])
print(data['vy'][1])
print(data['vz'][1])

