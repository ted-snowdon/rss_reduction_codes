import numpy as np
from matplotlib import pyplot as plt

from astropy import units as u
from astropy import modeling as mod
from astropy.io import fits
from astropy.table import Table

from specutils import Spectrum1D
from specidentify import WavelengthSolution
from specidentify import ArcIdentify

# enter the name of your arc spectrum here
#hdu = fits.open('cuar24125.fits')
#hdu = fits.open('cuar2525.fits')
hdu = fits.open('cuar26375.fits')
#hdu = fits.open('cuar32.fits')
#hdu = fits.open('cuar305.fits')
arc = Spectrum1D(spectral_axis=hdu[1].data['pixel']*u.pixel, flux=hdu[1].data['counts']*u.ct)

# a plot of the uncalibrated arc spectra
ax = plt.subplots()[1]  
ax.plot(arc.spectral_axis, arc.flux)  
ax.set_xlabel("Dispersion [pixels]")  
ax.set_ylabel("Flux")

# enter the name of your calibration line list file here
line_table = np.loadtxt('calib_cuar.txt', unpack=True, usecols=(0,1))
line_table = Table(line_table.T, names=('pixel', 'flux'))

xarr = np.arange(len(arc.data))

ws_init = mod.models.Chebyshev1D(3)
ws_init.domain = [xarr.min(), xarr.max()]
ws_init = mod.fitting.LinearLSQFitter()(ws_init, xarr, xarr)

ws = WavelengthSolution(x=None, wavelength=None, model=ws_init)

# w should be selected lines from the calibration file, x should be their locations in the arc spectrum
#24.125
#x = [607,673,732,789,1089,1541,1616,1756,1862,2191,2284,2356,2549,2668,2771,3064]
#w = [4072.385,4103.912,4131.724,4158.591,4300.101,4510.733,4545.052,4609.567,4657.901,4806.021,4847.81,4879.864,4965.08,5017.163,5062.037,5187.746]
#25.25
#x = [186,252,310,367,664,1113,1187,1326,1431,1665,1755,1847,1919,2227,2329,2616]
#w = [4072.385,4103.912,4131.724,4158.591,4300.101,4510.733,4545.052,4609.567,4657.901,4764.865,4806.021,4847.81,4879.864,5017.163,5062.037,5187.746]
#26.375
x = [238,684,757,895,1000,1231,1321,1413,1483,1672,1788,1888,2173,2890,3041]
w = [4300.101,4510.733,4545.052,4609.567,4657.901,4764.865,4806.021,4847.810,4879.864,4965.08,5017.163,5062.037,5187.746,5495.874,5524.957]
#32
#x = [72,160,238,313,713,1323,1425,1616,1762,2087,2215,2344,2445,2717,2886]
#w = [4072.385,4103.912,4131.724,4158.591,4300.101,4510.073,4545.052,4609.567,4657.901,4764.865,4806.021,4847.810,4879.864,4965.080,5017.163]
#30.5
#x = [631,720,799,875,1280,1898,2002,2197,2346,2678,2804,2942,3046]
#w = [4072.385,4103.912,4131.724,4158.591,4300.101,4510.733,4545.052,4609.567,4657.901,4764.865,4806.021,4847.810,4879.864]

ws = WavelengthSolution(x=x, wavelength=w, model=ws_init)
ws.fit()

# use wavelength 4100-5100 for 32, 3800-4900 for 30.5
# use 3800-5200 for 24.125, 4000-5400 for 25.25, 4200-5600 for 26.375
aw = ArcIdentify(arc, ws, line_table, wavelength_range=[4200, 5600], flatten_order=9)
aw.show_commands()
aw.draw_error()
aw.ws.sigma(aw.xp, aw.wp)

#hdu = fits.open('cuar24125.fits')
#hdu = fits.open('cuar2525.fits')
hdu = fits.open('cuar26375.fits')
#hdu = fits.open('cuar32.fits')
#hdu = fits.open('cuar305.fits')
agn = Spectrum1D(spectral_axis=aw.ws(hdu[1].data['pixel'])*u.angstrom, flux=hdu[1].data['counts']*u.ct)

# plot of the uncalibrated and calibrated arc spectra
ax = plt.subplots()[1]  
ax.plot(agn.spectral_axis, agn.flux, arc.spectral_axis, arc.flux)  
ax.set_xlabel("Dispersion [Angstrom]")  
ax.set_ylabel("Flux")

plt.show()

with open('calibrated_arc.dat', 'w') as file:
    for i in range(0, len(agn.spectral_axis)):
        wvl_item = str(agn.spectral_axis[i]).split(' ')[0]
        flx_item = str(agn.flux[i]).split(' ')[0]
        item_str = wvl_item + '\t' + flx_item + '\n'
        file.write(item_str)
file.close()
