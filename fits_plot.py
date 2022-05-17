import matplotlib.pyplot as plt
import os
import sys
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits

plt.style.use(astropy_mpl_style)

#example syntax: python fits_plot.py 20200101 10
# this will show mbxgpP2020010110.fits

date = sys.argv[1]
file = sys.argv[2]

folder = '/home/tedsnowdon/salt.data/'+date+'/product/'
filename = 'mbxgpP'+date+file+'.fits'

os.system('cp '+folder+filename+' /home/tedsnowdon/salt.reductions/titus_saures_rex/')

image_file = get_pkg_data_filename(filename)
image_data = fits.getdata(image_file, ext=1)

os.system('rm '+filename)

plt.figure()
plt.imshow(image_data, aspect='auto', cmap='inferno')
plt.clim(0.,100.)
plt.show()