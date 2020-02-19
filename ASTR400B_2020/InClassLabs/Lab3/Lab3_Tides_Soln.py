
# coding: utf-8

# In[1]:

# In Class Lab 3
# G. Besla 

# import relevant modules 
import astropy.units as u
import numpy as np


# The Large Magellanic Cloud is at a distance of 50 kpc from the Galactic Center. 
# It is observed to have a stellar disk that extends to a radius of at least 18.5 kpc.
# 
# ![LMC](./Lab3_Tidal.png)
# Deep photometric imaging reveals the faint stellar outskirts of the LMC. 
# Outskirts: DECam data Mackey+2016 MNRAS 459, 239. 
# Inner: shallower imaging from robotic telescopes Besla+2016 APJ 825.
# 
# In this lab we will determine
# the minimum mass required for the LMC so that it maintains the observed radius 
# in the face of the Milky Way's tidal field. 

# # Part A
# 
# We define the mass profile of the Milky Way using a Hernquist profile.
# 
# 
# $\rho(r) =  \frac{M_{halo}}{2\pi} \frac{a}{r(r+a)^3} \qquad M(r) =  \frac{M_{halo} r^2}{(a+r)^2}$ 
# 
# 

# ## #1
# 
# Create a function `HernquistM` that returns the dark matter halo mass at a given radius in units of solar mass.
# This function should take as input:  the distance from the Galactic center $r$, the scale radius $a$, and the halo mass $M_{halo}$.
# 
# 
# For the Hernquist scale radius for the Milky Way, use $a=60$ kpc. 
# 
# For $M_{halo}$ use your answer for the total mass of the simulated Milky Way you computed in Homework 3. 

# In[2]:

# Question 1
# Function that returns the Hernquist 1990 Mass profile
def HernquistM(r,a=60, Mhalo=1.96):
# Input:  Radius (kpc), Hernquist Scale Length (kpc), Total Mass (1e12 Msun)
# Defaults for scale length is 60 kpc, Default for Halo Mass is for MW from Assignment 3 
# Returns: Mass in units of Msun
    return round(Mhalo*r**2/(a+r)**2,2)*1e12*u.Msun


# In[65]:

print(HernquistM(10000))


# In[66]:

print(HernquistM(260))


# In[67]:

print(HernquistM(50))


# ## #2
# 
# Compute the total mass of the Milky Way within 50 kpc, including it's baryonic component (Dark Matter + Bulge + Disk), in units of solar mass. Use your answers from Homework 3 for the Bulge and Disk Masses. 
# Store this as a variable called `MassMW50`.
# 

# In[4]:

# MW disk mass and bulge mass
# using answers from assignment 3 
Mdisk =  7.5e10*u.Msun
Mbulge = 1e10*u.Msun


# In[5]:

# determine the total mass of the MW within 50 kpc
# Kochanek+1996 find 4.9e11
MassMW50 = HernquistM(50,60,1.96) + Mdisk + Mbulge
MassMW50


# # Part B
# 
# The Jacobi Radius for a satellite on a circular orbit about an extended host, where 
# the host is assumed to be well modeled as an isothermal sphere halo:
# 
# 
# $R_j = r  \bigg( \frac{M_{sat}}{2 M_{host}(<r)} \bigg)^{1/3}$
# 
# 
# The Isothermal Sphere approximation is not a bad one within 50 kpc.
# 
# Note also that the LMC is not on a circular orbit, but it is very close to its pericentric approach, where the velocity is all in the tangential component. So this isn't a terrible approximation either. 
# 
# ## #1
# Create a function called `JacobiMass` that returns the total mass of a satellite galaxy in units of Msun, 
# such that it has a given size 
# 
# Do this by rearranging the Jacobi Radius equation to solve for the satellite mass. 
# 

# In[52]:

# Part B . Question 1 

# Function that returns the minimum Satellite Mass
# needed to maintain a satellite's size and properties of a host galaxy.
# Given the Equation for the Jacobi Radius for an extended host, 
# we can rearrange the equation so that it returns the 
# minimum required satellite mass to maintain an observed radius
def JacobiMass(Rj,r,Mhost):
# input :  Rj = Jacobi Radius is approximated by the Observed Size of the Satellite in kpc
#        : r = distance to host in kpc
#        : Mhost = Host Mass enclosed within r in units of Msun
#  Returns: minimum satellite mass within its current size, in Msun
    return (Rj/r)**3*2*Mhost
    


# ## #2 
# 
# Determine the minimum total mass of the LMC needed to maintain its radius of 18.5 kpc in the face of the Milky Way's tidal 
# field at its current distance of 50 kpc. Store this as a variable called `LMCJacobiMass`.

# In[7]:

# Observed size of the LMC disk
# Mackey+2016
SizeL = 18.5*u.kpc
DistL = 50.0*u.kpc


# In[54]:

# LMC minimum mass in maximal MW halo model (from Simulation)
LMCJacobiMass = JacobiMass(SizeL,DistL,MassMW50)
print(LMCJacobiMass)


# ## #3
# 
# Recall that, ignoring centrifugal forces and assuming the host is a point mass, the tidal radius is given as :
# 
# $r_{tide} = r\left (\frac{m_{sat}}{4M_{host} } \right)^{1/3} $
# 
# Create a function to compute the total mass the must LMC possess to have a disk with radius 18.5 kpc.
# 
# The function should take as input, the current size of the satellite (kpc), this distnce ot the host(kpc) and the mass of the host (in Msun)
# 
# Use the function to determine the needed LMC mass and store it as a variable called `LMCTidalMass`. 

# In[19]:

# Function to compute the mass of the satellite needed such that the tidal radius is the current size of the 
# galaxy, ignoring centrifugal forces
def TidalMass(rtide, r, Mhost):
    # Input:  rtide, the tidal radius in kpc
    #        r, distance from host in kpc
    #       Mhost: Mass of the host galaxy in Msun
    #Returns: the satellite mass in Msun
    return 4*Mhost*(rtide/r)**3 


# In[21]:

LMCTidalMass = TidalMass(SizeL,DistL,MassMW50)
print(LMCTidalMass)


# ## #4
# 
# Compare `LMCTidalMass` to the calculation using the Jacobi Radius.
# How does the total mass of the LMC compare to its stellar mass (M$_\ast = 3 \times 10^9$ M$_\odot$)? 
# 

# In[25]:

print(LMCTidalMass/LMCJacobiMass)
# Because of centrifugal forces the minimum mass is smaller using the Jacobi Radius. 


# In[23]:

LMCMstar = 3e9*u.Msun


# In[26]:

print(LMCJacobiMass/LMCMstar)
# Mass ratio is ~ factor of 20. 


# # Part C: Consistency Check
# 
# 
# The maximal enclosed mass of the LMC at any radius can be determined by assuming a flat rotation curve. 
# 
# $V_c^2 = \frac{G M}{r} = constant$
#  
#  The rotation curve of the LMC is observed to flatten at a value of 91.7 +/- 18.8 km/s  (van der Marel & Kallivayalil 2014 ApJ 781)
#   (note that 1 kpc/Gyr $\sim$ 1 km/s). 
#  
#  Create a function called `MaxMass` that takes as input Vc (km/s) and distance to from the center of the galaxy (r) and returns the maximal dynamical mass in Msun. 
#  
# 

# In[27]:

# gravitational constant in units of kpc^3/Gyr^2/Msun
G = 4.498768e-6*u.kpc**3/u.Gyr**2/u.Msun


# In[45]:

# Maximal LMC mass
# Assuming LMC has a flat rotation curve to 18.5 kpc
# Vc = 91.7 +/- 18.8 km/s  van der Marel & Kallivayalil 2014
# 1 km/s ~ 1 kpc/Gyr
Vc = 91.7*u.kpc/u.Gyr


# In[37]:

# MLMC Vc^2 = GM/R = constant, rearrange for M:
def MaxMass(Vcirc=Vc, R=SizeL):
    return Vcirc**2/G*R


#  
# ##  #1  
# Compute the maximal mass enclosed by the LMC within the observed radius. Store it as a variable called `LMCMassMax`. 

# In[40]:

LMCMassMax = MaxMass()
print(LMCMassMax)


# ## #2
# 
# Is `LMCMassMax` consistent with `LMCJacobiMin`, the minimum mass needed to explain the observed radius of the LMC given the tidal field of the MW? If not, how could the numbers be reconciled?

# In[41]:

print(LMCMassMax/LMCJacobiMass)


# The minimum mass needed seems larger than the maximal mass possible.
# Either the LMC rotation curve needs to be higher (which it could within the errors)
# Or MW halo mass within 50 kpc is smaller, e.g. $3\times 10^{11}$ M$_\odot$.
# 
# Although note, the values are pretty close to being the same. 

# In[57]:

# Try increasing the LMC circular speed by 1 sigma
Vnew = Vc + 18.8*u.kpc/u.Gyr
print(MaxMass(Vcirc=Vnew))



# In[55]:

# LMC minimum mass in lower mass MW halo model  (recall Hernquist model gives 4e11 Msun within 50 kpc)
MinMW = 3e11*u.Msun
print(JacobiMass(SizeL,DistL,MinMW))


# Note that another equation can be used to describe the Jacobi radius if both the LMC and MW have flat rotation curves to large distances.

# In[ ]:

# From van der Marel et al. 2002 
# assuming a flat rotation curve for the MW and 
# for the satellite
def JacobiFlatVc(rsep,VcSat,VcHost):
    return rsep*(VcSat/VcHost)


# In[ ]:

# MW Vc = 206 to get M(<50 )= 4.9e11, assuming isothermal sphere.
JacobiFlatVc(DistL,91.7,240)


# In[ ]:



