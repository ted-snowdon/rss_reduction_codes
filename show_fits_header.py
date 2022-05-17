from astropy.io import fits

filename = input('Please enter fits filename: ')

hdul = fits.open(filename)

hdr = hdul[0].header

print('+---------- FITS INFO ----------+')
hdul.info()
print('+------------ HEADER -----------+')
print(repr(hdr))