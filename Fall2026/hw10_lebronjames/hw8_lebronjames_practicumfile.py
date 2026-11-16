#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Shiloh
#AST4762C
#Practicum 5/Homework 8


# In[2]:


"""
Researchers studied HD 189733b, and exoplanet referred to as a hot jupiter, to study how efficiently the tidally locked planet distributes heat from 
the perpetual dayside to the dark side. By observing the system's change in brightness as different sides of the planet rotate into view, the research
team was able to generate a longitudinal temperature map of the planet. The map showed a small temperature difference between the two sides of the
planet, unlike what is expected from hot jupiters. Interestingly, the brightest and hottest region of the planet was not at the point closest to the
star, rather, offset by nearly 30 degrees. From over 278,000 images, it was theorized that the planet has strong winds and efficient heat transport 
throughout, leading to the relatively uniform temperature observed.
"""


# In[3]:


"""
Researchers utilized the Spitzer Space Telescope to measure infrared emission from HD 189733b during a secondary eclipse. Results showed 
deviations from a normal blackbody curve, suggesting a low opacity region of the atmosphere, as well as traces of water vapor. It also demonstrates
that the planet does not have a temperature inversion in its upper atmosphere. Comparing the recorded spectrum to existing models, researchers 
hypothesize that HD 189773b does not efficiently redistribute heat energy to the dark side. This, however, could vary as the atmospheric
depth or wavelength of radiation changes, as suggested by tension in phase curve observations around 8 microns. These findings demonstrate
the complexity of the atmospheres of hot jupiters, and show how multi-wavelength eclipse observations can reveal a great deal of information
"""


# In[4]:


"""
Researchers observed HD 189733b using the Hubble Space Telescope's Spectroscopic instruments to measure how much light the planet reflects during
eclipses. Since thermal emission is negligible within the selected wavelength range, the reflected radiation is almost entirely starlight. This
allowed the team to directly determine the planet's geometric albedo. They found that the planet reflects substantial amounts of blue light, but 
decreases quickly as the wavelength increases. This phenomenon is best explained by the presence of reflective clouds or mist on the bright side
of the planet, as well as sodium absorption of reflected light beyond 450 nm. These two characteristics give the planet a lovely deep blue 
appearance, unlike anything we observe in our solar system
"""


# In[5]:


"""
Utilizing the James Webb Space Telescope's Mid-Infrared INstrument, researchers observed a complete infrared phase curve of HD 189733b within the
range 5-12 microns, allowing much more precise measurements than previously collected data. The resulting phase curve showed a much smaller
contrast between the temperatures of the two sides of the planet, showing clear signs of water vapor and powerful winds circulating heat across
the planet's surface. These observations reveal that the dark side of the planet is warmer than expected, and most likely sees high altitude silicate
clouds. Most importantly, the demonstrated a need a need to update older thermal models to account for the warm and cloudy atmosphere discovered
with the JWST, and showed just how impactful the telescope can be on exoplanet studies.
"""


# In[6]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = np.loadtxt("hw7_synthetic_transit_lightcurve.csv", delimiter=',', skiprows=1)

time = data[:,0]
flux = data[:,1]


# In[7]:


plt.figure(figsize=(8,5))
plt.scatter(time, flux, s=8, color="black")
plt.xlabel("Time (Days)", fontsize=12) #FIXED 18-Oct
plt.ylabel("Flux", fontsize=12) #FIXED 18-Oct
plt.title("Synthetic Transit Lightcurve")
#FIXED 18-Oct
plt.grid(alpha=0.24)
plt.tick_params(direction='in', labelsize=12)
plt.show()


# In[8]:


def roddyrich(x, center, width, depth):
    y=np.ones(len(x))
    zone = np.abs(x  - center) < (0.5 * width)
    y[zone] = 1-depth

    return y

guesscenter = np.median(time)
guesswidth = 0.02
guessdepth = 0.01

popt, pcov = curve_fit(roddyrich, time, flux, p0=[guesscenter, guesswidth, guessdepth])

supermodel = roddyrich(time, *popt) #FIXED 18-Oct
centerfit, widthfit, depthfit = popt
ucenter, uwidth, udepth = np.sqrt(np.diag(pcov))

print("Best center:", centerfit)
print("Best width:", widthfit)
print("Best depth:", depthfit)
print("Uncertainties:", ucenter, uwidth, udepth)
plt.scatter(time, flux, label='Observations')
plt.plot(time, supermodel, color='red', label='Best Fit Model') #FIXED
plt.grid(alpha=0.24)
plt.tick_params(direction='in', labelsize=12)
plt.xlabel("Time (Days)", fontsize=12) #FIXED 18-Oct
plt.ylabel("Flux", fontsize=12) #FIXED 18-Oct
plt.title("Synthetic Transit Lightcurve")
plt.legend()
plt.show()


# In[9]:


#curve_fit can't estimate covariance, function is essentially discontinuous


# In[10]:


radstar = 1.14
radplan = radstar * np.sqrt(depthfit)
print("Planet radius:", radplan)

d = depthfit
sigstar = 0.01
sigd = 0.002

uno = (np.sqrt(d) * sigstar)**2
dos = ((radstar / (2*np.sqrt(d))) * sigd)**2

sigplan = np.sqrt ( uno + dos ) 
print("Planetary Radius Uncertainty: ", sigplan)


# In[11]:


def normmedcomb(data, region=None):
    backup = data.copy()
    if region is None:
        normfactor = np.median(data, axis=(1,2))

    else:
        (y1, x1), (y2, x2) = region
        selected = data[:, y1:y2, x1:x2]
        normfactor = np.median(selected, axis=(1,2))
    newdata = backup / normfactor [:, None, None]
    combined = np.median (newdata, axis=0)
    return combined, normfactor


# In[12]:


from astropy.io import fits

hw6data = fits.getdata("darkdataforprac5.fits")
region = ((225,225), (-225, -225))
result = normmedcomb(hw6data, region)

combined= result[0]
normfactor = result[1]


# In[13]:


print("Normalization factors:", normfactor)
print("Shape of Combined Image:", np.shape(combined))
fits.writeto(
    "sky_13s_mednorm.fits",
    combined,
    overwrite=True,
)


# In[14]:


#Normalization factors: [10864.  11057.  10536.  10558.5 11046.  10623.   9916.5  9916.5  9859.5]
#Shape of Combined Image: (1024, 1024)


# In[15]:


#Shiloh
#homework 8
#AST4762C

def skycormednorm(object, sky, region=None):
    backobj = object.copy()
    backsky = sky.copy()
    if region is None:
        normfact = np.median(backobj)
    else:
        (y1,x1), (y2,x2) = region
        selected = backobj[y1:y2, x1:x2]
        normfact = np.median(selected)
    corrected = backobj - (backsky * normfact)
    
    
    return corrected


# In[19]:


region = ((225, 225), (-225, -225))

for i in range(9):
    corr = skycormednorm(hw6data[i],combined, region)

print(corr[217,184])


# In[20]:


fits.writeto(
    "stars_13s_9_nosky.fits",
    corr,
    overwrite=True,
)


# In[21]:


test = fits.getdata("stars_13s_9_nosky.fits")

print(np.shape(test))


# In[ ]:




