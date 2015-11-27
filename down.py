#downgrader.py

import numpy as np
from os import remove
import os.path
import astropy.io.fits as pyfits
import astropy
import astropy.constants
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline
from david_instrument import *
import pyfits
import matplotlib.pyplot as pyplot


header=pyfits.open('/Users/Daniel/Downloads/ic0719.fits')
flux=header[1].data['flux']
wavelength=header[1].data['loglam']
ivar = header[1].data['ivar']
error=np.sqrt(1/ivar)


r_instrument = np.zeros(len(wavelength))
for wi,w in enumerate(wavelength):
		if w<6000:
			r_instrument[wi] = (2270.0-1560.0)/(6000.0-3700.0)*w + 420.0 
		else:
			r_instrument[wi] = (2650.0-1850.0)/(9000.0-6000.0)*w + 250.0

r_instrument = np.zeros(len(wavelength))			
sres=r_instrument.fill(2.54)
log_wave=True
wave=wavelength

new_sres=sres/32.0
				
new_flux, matched_sres, sigma_offset, new_mask = \
        match_spectral_resolution(wave, flux, sres, wave, new_sres, min_sig_pix=0.0, log10=log_wave, new_log10=log_wave)
        
pyplot.plot(wave, new_flux, 'r')
pyplot.plot(wave, flux,'k')
pyplot.show()

