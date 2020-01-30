#!/usr/bin/env python
# coding: utf-8

# This code will read a file MW_000.txt. It will print out the time, total number of particles as variables, particle type, mass, x,y,z, and vx,vy,vz.

# In[22]:


#import libraries
#use abreviations to simplify later code
import numpy as np
import astropy.units as u


# In[23]:


filename = "MW_000.txt" #defines the file we will use

#define a function where we will read, open, and print data from the file MW_000.txt
def Read(filename):
    #Inputs:
    #filename, in this case MW_000.txt
    
    #Returns:
    #time, tot_part, data arrays from the file
    
    
    #opens and reads file
    file = open(filename,'r') 
    
    #reading line 1,labeling the value that is extracted, and storing the time
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*10.0*u.Myr
    
    #reading line 2,labeling the value that is extracted, and storing the total number of particles (tot_part)
    line2 = file.readline()
    label, value = line2.split()
    tot_part = float(value)
    
    #close the file
    file.close
    
    #store rest of file in np.genfromtxt to use the header information
    #parameters:
        #"dtype=None" splits line using white spaces
        #"names=True" creates arrays to store data with the right labels
            #labels are: type, m, x, y, z, vx, vy, vz
        #"skip_header=3" skips the first 3 lines
    data = np.genfromtxt("MW_000.txt", dtype=None, names=True, skip_header=3)
    
    return time, tot_part, data

#print the information form the file
if __name__=="__main__":
    
    #looks at the function above and extracts the data for time, tot_part, and data
    time,tot_part,data = Read(filename)
    
    #general print statements to print time, tot_part, and data
    print("the time is ", time)
    print("the total number of particles is ", tot_part)
    print(data)
    
    #selects piece of data to print
    print(data['type'[1]])
   


# In[ ]:




