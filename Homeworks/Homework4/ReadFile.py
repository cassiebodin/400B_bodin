#!/usr/bin/env python
# coding: utf-8

# This code will read a file MW_000.txt. It will print out the time, total number of particles as variables, particle type, mass, x,y,z, and vx,vy,vz.

# In[3]:


#import libraries
#use abreviations to simplify later code
import numpy as np
import astropy.units as u


# In[4]:


#define a function where we will read, open, and retrieve data from a file designated as "filename"
def Read(filename):
    #Inputs:
    #filename, we will use MW_000.txt to test the code
    
    #Returns:
    #time (Myrs)
    #tot_part, total number of particles ()
    #data arrays from the file
    
    
    #opens and reads file
    file = open(filename,'r') 
    
    #reading line 1,storing the time
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr
    
    #reading line 2, storing the total number of particles (tot_part)
    line2 = file.readline()
    label, value = line2.split()
    tot_part = float(value)
    
    #close the file
    file.close
    
    #store rest of file in np.genfromtxt to use the header information
    #parameters:
        #"dtype=None" specifies data type, none is default float
        #    default deliminater splits line using white spaces
        #"names=True" creates arrays to store data with the right labels
            #labels are: type, m, x, y, z, vx, vy, vz
        #"skip_header=3" skips the first 3 lines
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)
    
    #returns the time, total number of particles, and array of data
    return time, tot_part, data


# In[5]:


#test code for the file MW_000.txt
time, tot_part, data = Read("MW_000.txt")


# In[6]:


#print time
time


# In[7]:


#print total number of particles
tot_part


# In[8]:


#prints mass of the first particle, note: first particle is particle 0
data['m'][0]*u.Msun*1e10


# In[9]:


#print type of first particle, note: first particle is particle 0
data['type'][0]


# In[ ]:




