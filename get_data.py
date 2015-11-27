#get_data.py

import numpy as np
import pandas as pd
import astropy.io.fits as pyfits
from IPython.core.debugger import Tracer
from dust import *
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord

def get_data(options):

	"""
	Retrieves data from the location and format specified,
	returns data SED, error SED, metadata.
	"""

	if options['data_type']=='custom':

		if options['data_format'] == 'ascii':
			"""
			Format of ASCII files MUST be:
			0 lines of header
			3 columns of numeric wavelength, flux and error.
			For custom data formats, the following information is also provided in
			the parameters file: redshift, vdisp,ra,dec (if not required, set to -99)
			"""

			#import custom_properties
			wavelength,flux,error   = np.loadtxt(options['file_in'], skiprows=0, unpack=True)
			bad_flags               = np.ones(len(wavelength))
			ra          = options['ra']
			dec         = options['dec']
			redshift    = options['redshift']
			vdisp       = options['vdisp']  
			trust_flag 	= 1  
			r_instrument = wavelength/3.0 
			objid = ''       

		elif options['data_format'] == 'ascii_claire':
			"""
			Special ASCII files with redshift in the file name in this format:
			(for collaboration work with Claire Le Cras)

			0 lines of header
			3 columns of numeric wavelength, flux and error.
			For custom data formats, the following information is also provided in
			the parameters file: redshift, vdisp,ra,dec (if not required, set to -99)
			"""

			#import custom_properties
			wavelength,flux,error   = np.loadtxt(options['file_in'], skiprows=0, unpack=True)

			wavelength = wavelength*(1.0+redshift)
			bad_flags               = np.ones(len(wavelength))

			redshift = (float(options['file_in'].split('/')[-1].split('_')[4])+\
						float(options['file_in'].split('/')[-1].split('_')[6]))/2.0
			ra          = 0#options['ra']
			dec         = 0#options['dec']
			redshift    = options['redshift']
			vdisp       = 170#options['vdisp']  
			trust_flag 	= 1      

			wavelength = wavelength*(1.0+redshift)

		if options['data_format'] == 'fits_matt':
			"""
			Format of fits files MUST be:
			- wavelength encoded as SDSS in linear units, e.g:
			CDELT1  =     2.10055088996887
			CRVAL1  =     3406.70236811121
			CRPIX1  =    -162.432189873714

			A file for spectra, and a file for error,
			with ***s.fits and ***e.fits as file names.
			"""
			#NGC_891_test.ms.fits

			hdus = pyfits.open(options['file_in'])
			crval1      = hdus[0].header['CRVAL1']
			cdelt       = hdus[0].header['CDELT1']
			naxis       = hdus[0].header['NAXIS1']
			crpix=hdus[0].header['CRPIX1']

			ra_str      = hdus[0].header['TARGRA']
			dec_str     = hdus[0].header['TARGDEC']

			radec = SkyCoord(ra_str+' '+dec_str, unit=(u.hourangle, u.deg))
			ra,dec = radec.ra.value,radec.dec.value
			redshift    = 0.0001761#hdus[0].header['Z_FIN']
			crval=crval1-(crpix*cdelt)
			wavelength=np.arange(naxis)*cdelt+crval

			
			flux  = hdus[0].data[options['file_hdu']]*10.0**(-17)
			ehdu = pyfits.open(options['file_in'][:-6]+'e.fits')
			error  = ehdu[0].data[options['file_hdu']]*10.0**(-17)

			bad_flags               = np.ones(len(wavelength))
			vdisp = options['vdisp']
			trust_flag=1
			objid = hdus[0].header['OBJECT']
			r_instrument = np.zeros(len(wavelength))
			r_instrument[:] = 2000

			pass

	elif options['data_type'] == 'sdss_dr7':
		hdus = pyfits.open(options['file_in'])

		crval1      = hdus[0].header['CRVAL1']
		cdelt       = hdus[0].header['CD1_1']
		naxis       = hdus[0].header['NAXIS1']
		crpix=0.
		redshift    = hdus[0].header['Z_FIN']
		crval=crval1-(crpix*cdelt)
		initial_initial_data_wavelength=10**(np.arange(naxis)*cdelt+crval)

		# In SDSS FITS data, generally HDU0 = flux, HDU5 = error, HDU4 = good pixel flag,
		# HDU=3 is e-cleaned 
		flux  = hdus[0].data#*hdus[4].data
		error = hdus[5].data#*hdus[4].data

		bad_flags   = hdus[4].data
		#mean_signal = np.sum(initial_data_flux)/np.size(initial_data_flux)

		wavelength  = initial_initial_data_wavelength
 		header_test=pyfits.open('/Users/Daniel/firefly_public/DR7_sigma_catalogue.fits')
		main_header=header_test[1].data

		sdss_plate=[]
		mjd=[]
		sdss_fiber=[]
		velocity_dispersion=[]
		ra_1=[]
		dec_1=[]
		z=[]
		
		for i in range(len(main_header)):
			element=main_header[i]
		 	plate=element[0]
		 	mj=element[1]
		 	fiber=element[2]
		 	sigma=element[3]
			right_asc=element[4]
			declination=element[5]
		 	redshi=element[6]
		 	sdss_plate.append(plate)
		 	mjd.append(mj)
		 	sdss_fiber.append(fiber)
		 	velocity_dispersion.append(sigma)
		 	z.append(redshi)
			ra_1.append(right_asc)
			dec_1.append(declination)
	
		sdss_plate=np.array(sdss_plate)
		mjd=np.array(mjd)
		sdss_fiber=np.array(sdss_fiber)
		velocity_dispersion=np.array(velocity_dispersion)
		ra_=np.array(ra_1)
		dec_1=np.array(dec_1)
		z=np.array(z)
		 
		trial=options['file_in'].split('-')
		trial_plate=int(trial[1])
		trial_mjd=int(trial[2])
		trial2=trial[3].split('_')
		trial_fiber=int(trial2[0])
		ra=0.0
		dec=0.0

		
		index=np.where((sdss_plate == trial_plate) & (mjd == trial_mjd) & (sdss_fiber == trial_fiber))
		if not velocity_dispersion[index]:
			vdisp=0
		else:
		 	vdisp=float(velocity_dispersion[index])
		if velocity_dispersion[index] == -9999.0:
			vdisp=0
		
		ra=0.0
		dec=0.0
 		trust_flag=1
		objid = options['file_in']
 		 
 		r_instrument = np.zeros(len(wavelength))
		for wi,w in enumerate(wavelength):
			if w<6000:
				r_instrument[wi] = (2270.0-1560.0)/(6000.0-3700.0)*w + 420.0 
			else:
				r_instrument[wi] = (2650.0-1850.0)/(9000.0-6000.0)*w + 250.0


	elif options['data_type'] == 'sdss_boss':
		pass
	elif options['data_type'] == 'm67':
		# M67-specific stuff:
		m67_file = 'data/m67/int_m67_x135.txt'

		data = open(m67_file, 'r')
		initial_data_lambda_list,initial_data_flux_list = [],[]

		for i,line in enumerate(data.readlines()):
			initial_data_lambda_element,initial_data_flux_element = line.split()
			initial_data_lambda_list.append(float(initial_data_lambda_element))
			initial_data_flux_list.append(float(initial_data_flux_element))

		wavelength = np.asarray(initial_data_lambda_list[1:1000]+initial_data_lambda_list[999:])
		flux 		= np.asarray(initial_data_flux_list) #/ initial_data_wavelength
		redshift 			= 0.0001
		trust_flag 			= 1
		vdisp			= 5.0

		bad_flags = np.zeros(len(flux))
		bad_flags = bad_flags+1
		error 		= flux /100.0
		ra=0.0
		dec=0.0
		error[np.where(flux==0)]=np.median(flux)*10000000.0
		bad_flags[np.where(flux==0)] = 0
		objid = options['file_in']
		r_instrument = np.zeros(len(wavelength))
		r_instrument[:] = 10000


	elif options['data_type'] == 'ngc':
		pass
	elif options['data_type'] == 'mock_ssp':
		pass
	elif options['data_type'] == 'mock_csp':
		pass
	elif options['data_type'] == 'moresco_boss':

		file_split = options['file_in'].split('_')
		print file_split
		vdisp_str = file_split[2]
		if vdisp_str == 'mm1':
			vdisp = 100.0#*1.5
		elif vdisp_str == 'mm2':
			vdisp = 175.0#*1.5
		elif vdisp_str == 'mm3':
			vdisp = 225.0#*1.5
		elif vdisp_str == 'mm4':
			vdisp = 275.0#*1.5
		elif vdisp_str == 'mm5':
			vdisp = 400.0
		print vdisp
		redshift = float(file_split[4][1:])
		wavelength,flux,error_part1,junk1,junk2,error_part2   = \
						np.loadtxt(options['file_in'], skiprows=1, unpack=True)

		error = error_part1 / np.sqrt(error_part2)
		wavelength = wavelength*(1.0+redshift)
		ra=0.0
		dec=0.0
		bad_flags = np.ones(len(flux))
		error[np.where(flux==0)]=np.median(flux)*10000000.0
		bad_flags[np.where(flux==0)] = 0
		trust_flag = 1
		objid = options['file_in']



		r_instrument = np.zeros(len(wavelength))
		for wi,w in enumerate(wavelength):
			if w<6000:
				r_instrument[wi] = (2270.0-1560.0)/(6000.0-3700.0)*w + 420.0
			else:
				r_instrument[wi] = (2650.0-1850.0)/(9000.0-6000.0)*w + 250.0

		# grad = (2500.0-1500.0)/(9000-3800)
		# r_instrument = wavelength*grad + 770

		# kyle_file = np.zeros((3,len(wavelength)))
		# kyle_file[0] = wavelength/(1.0+redshift)
		# kyle_file[1] = flux
		# kyle_file[2] = r_instrument
		# np.savetxt('vdisp_test_for_kyle.txt',kyle_file.T)


		# error = error + np.median(flux)*0.1
		# plt.plot(wavelength,flux,'k')
		# plt.plot(wavelength,error,'cyan')
		# plt.show()


	elif options['data_type'] == 'test_manga' or options['data_type'] == 'manga':


		header      = pyfits.open(options['file_in'],ignore_missing_end=True)
		flux_all    = header['FLUX'].data
		r_instrument = header['SRES'].data
		emission    = header['ELOMEW'].data
		wavelength  = header['WAVE'].data
		ivar        = header['IVAR'].data    
		length      = np.shape(flux_all)[0]
		maxshape    = np.shape(flux_all)[1]

		bin_number  = options['bin_number']

		flux_em     = flux_all[:,bin_number]
		em_line     = emission[:,bin_number]
		flux        = flux_em-em_line
		inv_variance= ivar[:,bin_number]
		error       = np.sqrt(1.0/(inv_variance))

		ex_data = header[0].header
		ra      = ex_data['OBJRA']
		dec     = ex_data['OBJDEC']
		objid   = ex_data['MANGAID']
		pos_data= header[1].data
		xpos    = pos_data['XPOS'][bin_number]
		ypos    = pos_data['YPOS'][bin_number]


		#redshift = 0.0
		bad_flags = np.ones(len(flux))
		z_arr 	= header['STFIT'].data['KIN'][:,0] / 299792.458
		sig_arr = header['STFIT'].data['KIN'][:,1]
		vdisp = sig_arr[bin_number]
		redshift = z_arr[bin_number]
		ebv_mw = 0.0
		trust_flag = True

		print "REDSHIFT:"
		print redshift
		# MaNGA galaxies need to be re-redshifted and then de-redshifted
		# so that we still take account of luminosity distance later on.
		#wavelength = wavelength*(1.0+redshift)

	#else:
	#	print "data_type not recognised...stopping"
	#	raise NameError('Input data type not recognised... stopping. Edit parameter "data_type".')


	restframe_wavelength = wavelength / (1.0+redshift)

	if options['milky_way_reddening']:
		# Find amount of MW reddening on the models
		ebv_mw 				 	= get_dust_radec(ra,dec,'ebv')
  	else:
  		ebv_mw = 0.0

  	bad_data 			= np.isnan(flux) | np.isinf(flux) | (flux <= 0.0) | np.isnan(error) | np.isinf(error)
  	flux[bad_data] 		= 0.0
  	error[bad_data] 	= np.max(flux) * 99999999999.9
  	bad_flags[bad_data] = 0

  	if options['plot_diagnostics']:
		plt.plot(wavelength,flux)
		plt.show()
			
  	data_dict = {'wavelength':restframe_wavelength,'flux':flux,'error':error,\
			'flags':bad_flags,'redshift':redshift,\
			'vdisp':vdisp,'ra':ra,'dec':dec,'ebv_mw':ebv_mw,'trust_flag':trust_flag,'objid':objid,\
			'r_instrument':r_instrument}

	if options['observation_type']=='ifu':
		data_dict['bin_number'] = bin_number
		data_dict['xpos'] 	= xpos
		data_dict['ypos'] 	= ypos

	return data_dict





