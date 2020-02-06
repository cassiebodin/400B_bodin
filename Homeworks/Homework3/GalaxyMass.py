#!/usr/bin/env python
# coding: utf-8

# # Homework 3
# 
# In this code I will compute the mass breakdown of the Local Group (SnapNumber 0) using the most massive members: the Milky Way (MW), M31, and M33.

# ## 1 Return Mass
# 
# Create a function that will return the total mass of each of the galaxy components. 

# In[2]:


#import libraries
#use abreviations to simplify later code
import numpy as np
import astropy.units as u
from ReadFile import Read


# In[3]:


#define a function to calculate the total mass of the galaxy components

def ComponentMass(filename, par_type):
    #Inputs:
    #filename, designates which file we will take data from
    #    files we will use are: MW_000.txt, M31_000.txt, M33_000.txt
    #par_type, the particle type
    #     Halo:1
    #     Disk:2
    #     Bulge:3
    
    #Returns:
    #Mtot, the total mass of a galaxy component
    
    #pull data from the ReadFile using our define import
    time, tot_part, data = Read(filename)
    
    #creates an index for all the particles given a property
    index = np.where(data['type'] == par_type)
    
    #stores those components in the index as new variables 
    #mass in the file was labeled as in 10^10 Msun so multiply by 1e10
    mnew=data['m'][index]*1e10*u.Msun
    
    #generate the mass
    #Note: will round the mass to 3 decimal places later in the code
    Mtot= np.sum(mnew)
    
    #return the total mass
    return Mtot


# ## 2 Mass Breakdown
# 
# Make a table of the masses computed using the Component Mass function and the files for the Milky Way(MW) , M31, and M33. Files we will use are: MW_000.txt, M31_000.txt, M33_000.txt

# ### 2.1 Mtot for each type of particle (Halo, Disk, Bulge) in each Galaxy: MW, M31, and M33

# In[4]:


#Get data for the milky way file "MW_000.txt" for the total masses
#the different numbers after mtot identify the particle type as following:
#     Halo:1
#     Disk:2
#     Bulge:3
mtotmw_1 = np.round(ComponentMass("MW_000.txt",1)/1e12,3)
mtotmw_2 = np.round(ComponentMass("MW_000.txt",2)/1e12,3)
mtotmw_3 = np.round(ComponentMass("MW_000.txt",3)/1e12,3)

#Print the masses for each in 10^12 Msun
print("Total mass in the MW halo is ", mtotmw_1)
print("Total mass in the MW disk is ", mtotmw_2)
print("Total mass in the MW bulge is ", mtotmw_3)


# In[5]:


#Get data for the galaxy M31 file "M31_000.txt" for the total masses
#the different numbers after mtot identify the particle type as following:
#     Halo:1
#     Disk:2
#     Bulge:3
mtotm31_1 = np.round(ComponentMass("M31_000.txt",1)/1e12,3)
mtotm31_2 = np.round(ComponentMass("M31_000.txt",2)/1e12,3)
mtotm31_3 = np.round(ComponentMass("M31_000.txt",3)/1e12,3)

#Print the masses for each in 10^12 Msun
print("Total mass in the halo is ", mtotm31_1)
print("Total mass in the M31 disk is ", mtotm31_2)
print("Total mass in the M31 bulge is ", mtotm31_3)


# In[6]:


#Get data for the galaxy M33 file "M33_000.txt" for the total masses
#the different numbers after mtot identify the particle type as following:
#     Halo:1
#     Disk:2
#     Bulge:3
mtotm33_1 = np.round(ComponentMass("M33_000.txt",1)/1e12,3)
mtotm33_2 = np.round(ComponentMass("M33_000.txt",2)/1e12,3)
mtotm33_3 = np.round(ComponentMass("M33_000.txt",3)/1e12,3)

#Print the masses for each in 10^12 Msun
print("Total mass in the M33 halo is ", mtotm33_1)
print("Total mass in the M33 disk is ", mtotm33_2)
print("Total mass in the M33 bulge is ", mtotm33_3)


# ### 2.2 Mtot in each Galaxy: MW, M31, and M33

# In[7]:


#Calculate the total mass of the Milky Way by adding the mass components
MtotMW= mtotmw_1+mtotmw_2+mtotmw_3

#Print the mass in 10^12 Msun
print("Total mass in the MW is ", MtotMW)


# In[8]:


#Calculate the total mass of the M31 by adding the mass components
MtotM31= mtotm31_1+mtotm31_2+mtotm31_3

#Print the mass in 10^12 Msun
print("Total mass in M31 is ", MtotM31)


# In[9]:


#Calculate the total mass of the M33 by adding the mass components
MtotM33= mtotm33_1+mtotm33_2+mtotm33_3

#Print the mass in 10^12 Msun
print("Total mass in M33 is ", MtotM33)


# ### 2.3 Mtot of the local group

# In[10]:


#Calculate the total mass of the local group
MtotLG= MtotMW+MtotM31+MtotM33

#Print the mass in 10^12 Msun
print("Total mass in the local group is ", MtotLG)


# ### 2.4 Baryon Fraction for each Galaxy: MW, M31, and M33
# baryon fraction = total stellar mass/total mass = (mass in disk(type2)+mass in bulge(type3))/total mass
# or
# fbar=(mtot_2+mtot_3)/Mtot

# In[16]:


#Calculate the baryon fraction in the Milky Way
#need to calculate the stellar mass first
mstel_mw=np.round(mtotmw_2+mtotmw_3,3)

#then calculate baryon fraction
fbarMW=np.round((mstel_mw/MtotMW),3)

#print baryon fraction
print("baryon fraction for the milky way is", fbarMW)


# In[17]:


#Calculate the baryon fraction in the M31
#need to calculate the stellar mass first
mstel_m31=np.round(mtotm31_2+mtotm31_3,3)

#then calculate baryon fraction
fbarM31=np.round((mstel_m31/MtotM31),3)

#print baryon fraction
print("baryon fraction for M31 is", fbarM31)


# In[18]:


#Calculate the baryon fraction in the M33
#need to calculate the stellar mass first
mstel_m33=np.round(mtotm33_2+mtotm33_3,3)

#then calculate baryon fraction
fbarM33=np.round((mstel_m33/MtotM33),3)

#print baryon fraction
print("baryon fraction for M33 is", fbarM33)


# In[19]:


#for the homework print the stellar masses
print("stellar mass in mw",mstel_mw)
print("stellar mass in m31",mstel_m31)
print("stellar mass in m33",mstel_m33)


# ### 2.4 Baryon Fraction for the local group
# baryon fraction = total stellar mass/total mass = (mass in disk(type2)+mass in bulge(type3))/total mass
# or
# fbar=(mtotmw_2+mtot31_2+mtot33_2+mtotmw_3+mtotm31_3+mtotm33_3)/MtotLG

# In[15]:


#Calculate the baryon fraction in the local group
#found the total stellar mass first
MTotStel=mstel_mw+mstel_m31+mstel_m33

fbarLG=np.round((MTotStel/MtotLG),3)

#print baryon fraction
print("baryon fraction for the local group is", fbarLG)


# In[ ]:




