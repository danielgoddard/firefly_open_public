import numpy as np
from fitter import fitter
import seaborn as sns

"""

NOTE: This set of convenience routines are for developer use only!
For testing that FIREFLY is running correctly, run firefly in "test" mode (see README).

###

A variety of useful testing tools to make sure the FIREFLY fitter runs properly.

"""

def test_fitter(num_points):
	wavelength_in 	= np.arange(num_points)
	data_in 		= np.random.rand(num_points)*10.0
	error_in		= np.ones(num_points)*5.0

	models_in 		= np.zeros((100,num_points))

	for m in range(len(models_in)):
		models_in[m,:] = np.random.rand(num_points)*10.0


	fitter(wavelength_in,data_in,error_in,models_in)

def test_fitter_mocks():

	import pyfits
	import random
	import os

	metal_used 			= ['z001','z002','z004','z0001.bhb','z0001.rhb','z10m4.bhb','z10m4.rhb']
	age_used 			= [	'3M','3_5eM','4M','4_5eM','5M','5_5eM','6M','6_5eM',\
							'7M','7_5eM','8M','8_5eM','9M','9_5eM',\
							'10M','15M','20M','25M','30M','35M','40M','45M','50M','55M',\
							'60M','65M','70M','75M','80M','85M','90M','95M',\
							'100M','200M','300M','400M','500M','600M','700M','800M','900M',\
							'1G','1G_500M','2G','3G','4G','5G','6G','7G','8G','9G',\
							'10G','11G','12G','13G','14G','15G']


	sigma_use 			= '100'
	model_used_array 	= ['MILES_UVextended']
	signal_to_noise 	= 20.0

	initial_model_flux	= []
	save_norm_models	= []

	on = True
	for j in range(len(metal_used)):
		for i in range(len(age_used)):
			input_data = metal_used[j]+'/CONVERTED_'+age_used[i]+'_log'
			file_open = '/Users/david/Downloads/ss/'+input_data+'.fits'
			test_file = os.path.isfile(file_open)
			if test_file == False:
				print "CANNOT FIND FILE"
				print file_open
				trust_flag = 0
  				null_output_structure = {'flux':[0,0],'error':[0,0],'wavelength':[0,0],\
  										 'redshift':0,'sigma_in':0,'eline_correction':0,'trust_flag':0,\
  										 'out_correction':1}
  			else:
				hdus = pyfits.open(file_open)

				crval1 		= hdus[0].header['CRVAL1']
				cdelt  		= hdus[0].header['CD1_1']
				naxis  		= hdus[0].header['NAXIS1']
				crpix=0.	   
  				crval=crval1-(crpix*cdelt)

				initial_data_wavelength 	= 10**(np.arange(naxis)*cdelt+crval)
				flux_int 	= hdus[0].data

				
				initial_data_flux_list  = []
				initial_data_error_list = []
				bad_flags = np.zeros(len(initial_data_wavelength))
				bad_flags = bad_flags + 1

				save_norm_models.append(np.sum(flux_int[100:-20]))
				initial_model_flux.append(np.asarray(flux_int[100:-20])/np.sum(flux_int[100:-20]))

	
	

	# CSP test:
	file_open='/Users/david/Downloads/log_miles_kr_tau_1_dust_no_age_5.fits'
	print file_open
	
	test_file = os.path.isfile(file_open)
	if test_file == False:
		trust_flag = 0
		null_output_structure = {'flux':[0,0],'error':[0,0],'wavelength':[0,0],\
								 'redshift':0,'sigma_in':0,'eline_correction':0,'trust_flag':0,\
								 'out_correction':1}
	else:
		hdus 	= pyfits.open(file_open)
		data_in = hdus[1].data

		initial_data_wavelength = 10**np.asarray(data_in['wavelength'])[0]
		flux_int 		= np.asarray(data_in['flux'])

		signal_to_noise = 20.0 # typical BOSS x 10
		initial_data_flux_list  = []
		initial_data_error_list = []

		for i in range(len(flux_int[0])):
			random.seed()
			initial_data_flux_list.append(random.gauss(flux_int[0][i],flux_int[0][i]/signal_to_noise))
			initial_data_error_list.append(flux_int[0][i] / signal_to_noise)


		initial_data_flux = np.asarray(initial_data_flux_list)[100:-20]
		initial_data_error= np.asarray(initial_data_error_list)[100:-20]






	# MULTI-SSP test:
	# initial_data_flux 	= np.zeros(len(initial_model_flux[0]))
	# initial_data_error 	= np.zeros(len(initial_model_flux[0]))
	# for j in range(len(metal_used)):
	# 	for i in range(len(age_used)):
	# 		if metal_used[j] == 'z002' and age_used[i] in ['10G','1G','5G','7G','500M']:
	# 			input_data = metal_used[j]+'/CONVERTED_'+age_used[i]+'_log'
	# 			file_open = '/Users/david/Downloads/kr/'+input_data+'.fits'
	# 			test_file = os.path.isfile(file_open)
	# 			if test_file == False:
	# 				print "CANNOT FIND FILE"
	# 				print file_open
	# 				trust_flag = 0
	#   				null_output_structure = {'flux':[0,0],'error':[0,0],'wavelength':[0,0],\
	#   										 'redshift':0,'sigma_in':0,'eline_correction':0,'trust_flag':0,\
	#   										 'out_correction':1}
	#   			else:
	# 				hdus = pyfits.open(file_open)

	# 				crval1 		= hdus[0].header['CRVAL1']
	# 				cdelt  		= hdus[0].header['CD1_1']
	# 				naxis  		= hdus[0].header['NAXIS1']
	# 				crpix=0.	   
	#   				crval=crval1-(crpix*cdelt)

	# 				initial_data_wavelength 	= 10**(np.arange(naxis)*cdelt+crval)
	# 				flux_int 	= hdus[0].data

					
	# 				initial_data_flux_list  = []
	# 				initial_data_error_list = []
	# 				bad_flags = np.zeros(len(initial_data_wavelength))
	# 				bad_flags = bad_flags + 1


	# 				# for i in range(len(flux_int)):
	# 				# 	random.seed()
	# 				# 	initial_data_flux_list.append(random.gauss(flux_int[i],flux_int[i]/signal_to_noise))
	# 				# 	initial_data_error_list.append(flux_int[i]/signal_to_noise)
	# 				# 	if on:
	# 				# 		initial_data_flux = np.asarray(initial_data_flux_list)
	# 				# 		initial_data_error= np.asarray(initial_data_error_list)
	# 				# 		on = False
	# 				# 	else:
	# 				# 		print np.shape(np.asarray(initial_data_flux_list))
	# 				# 		print np.shape(initial_data_flux)

	# 				initial_data_flux = initial_data_flux + flux_int[100:-20]
					



				
	initial_data_error= initial_data_flux/signal_to_noise
	models = np.asarray(initial_model_flux)

	norm_models = np.asarray(save_norm_models)
	norm_factor = 1.0/(norm_models/np.sum(initial_data_flux))

	error = initial_data_error/np.sum(initial_data_flux)
	data = initial_data_flux/np.sum(initial_data_flux)
	
	# plt.plot(data)
	# plt.plot(error)
	# plt.show()

	weights,chis,branch = fitter(initial_data_wavelength[100:-20],data,error,models)

	best_sol 		= np.argmin(chis)
	mass_estimate 	= np.dot(norm_factor,weights[best_sol])
	all_masses 		= np.dot(norm_factor,weights.T)

	ind_sort = np.argsort(chis)
	# plt.plot(all_masses,chis,'o')
	# plt.show()
	# Tracer()()
	# Convert chis to probs:
	# assume normal distro with best fit = mean
	# hence variance is sqrt(2*mean)
	min_chis = np.min(chis)
	var_derive = 2*min_chis
	prob = np.exp((min_chis-chis)/2/var_derive)

	sns.tsplot(walks, ci=100, color=pal)

	plt.plot(all_masses,prob,'o')
	plt.show()

# import cProfile
# cProfile.run('test_fitter_mocks()')

# Tracer()()
time_begin = time.time()
test_fitter_mocks()
print "Time taken (seconds):"
print time.time()-time_begin
#test_fitter(10)