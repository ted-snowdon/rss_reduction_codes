# Plots graph of reduced fits data to check reduction was successful
# example syntax: python fits_check.py 20200101 10
# this will show 20200101/r010.fits

import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits

date = sys.argv[1]
frame = sys.argv[2].zfill(3)

filepath = '/home/tedsnowdon/salt.reduced/'+date+'/r'+frame+'.fits'

with fits.open(filepath) as hdul:
    data = hdul[0].data

plt.plot(data[0], data[1])
plt.plot(data[0], data[2])
plt.show()