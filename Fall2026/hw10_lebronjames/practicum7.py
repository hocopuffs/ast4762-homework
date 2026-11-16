#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np

def disk(r, cyx, image_2d):
    y, x = np.ogrid[:image_2d[0], :image_2d[1]]
    rad_mask = ((y - cyx[0])**2 + (x - cyx[1])**2) < r**2
    return rad_mask

def apphot(image, mask, portal, inner, outer, center):
    y0, x0 = center
    ny, nx = image.shape
    ymin = max(0, int(y0 - (outer + 5)))
    ymax = min(ny, int(y0 + (outer + 5) + 1))
    xmin = max(0, int(x0 - (outer + 5)))
    xmax = min(nx, int(x0 + (outer + 5) + 1))
    subimage = np.copy(image[ymin:ymax, xmin:xmax])
    submask = mask[ymin:ymax, xmin:xmax]

    yadj = y0 - ymin
    xadj = x0 - xmin

    outerdisk = disk(outer, (yadj, xadj), subimage.shape)
    innerdisk = disk(inner, (yadj, xadj), subimage.shape)
    skymask = outerdisk & ~innerdisk

    meanie = np.nanmean(subimage[skymask & (submask >0)])

    subimage -= meanie

    GlaDos = disk(portal, (yadj, xadj), subimage.shape)

    flux = np.nansum(subimage[GlaDos & (submask > 0)])

    #for bonus
    disgrace = np.sum(GlaDos & (submask ==0))

    return flux, meanie, disgrace


# In[8]:


size = 21
y, x = np.indices((size, size))
y0, x0 = 10, 10
sig = 2
gaussian = np.exp(-((x-x0)**2 + (y-y0)**2)/(2*sig**2))

mask = np.ones_like(gaussian)

flux, sky, disgrace = apphot(gaussian, mask, 6, 7, 10, (y0, x0))

print(flux, sky, disgrace)


# In[ ]:




