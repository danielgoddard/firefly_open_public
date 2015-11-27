import numpy as np
from IPython.core.debugger import Tracer
from parameters_obtain import *
from get_data import *
from get_model import *
from fitter import *
from normalise_spec import *
from scipy.stats import sigmaclip

import matplotlib.pyplot as plt
from match_data_models import *
from dust import *
from statistics import *
from astropy.cosmology import Planck13
import os.path

from star_formation_rate import star_formation_rate

def firefly_single(parameters):

	"""
	The routine for a single run of FIREFLY.
	It is called from firefly_job, test_firefly, or can be run 
	interactively for a custom single SED.
	
	In the interactive case, one sets the 'custom' and 'interactive' parameters
	in the parameter file, then one enters at an interative prompt:
	> from firefly_single import firefly_single
	> firefly_single('[locationpathofdata]','custom/[outputdirname]','./parameters.dat')
	One is then able to view the output plots at the native X window. 

	This routine retrieves the options from the parameters file,
	including the locations of the data file and models to be used.

	It then opens the data file, model files, matches their resolutions
	(downgrading to instrumental+velocity disp resolution if necessary)
	fits, then produces output files and plots. 

	INPUTS: 
	- options_file: location of the parameter file (default: ./parameters.dat) 

	No outputs.
	"""

    # data_file, output_dir
	

	data = get_data(parameters)

	# restrict_ages can be default (allow age<+1 Gyr age uni), off (allow all ages), strict (age<age uni only)
	if parameters['restrict_ages'] == 'default':
		age_universe = Planck13.age(data['redshift'])
		parameters['age_limits'][1] = np.log10(age_universe.value+1.0)+9.0 # log(yr units)
	elif parameters['restrict_ages'] == 'strict':
		age_universe = Planck13.age(data['redshift'])
		parameters['age_limits'][1] = np.log10(age_universe.value)+9.0

	# Get the models with observation information needed for downgrading.
	for mi,mm in enumerate(parameters['model_libs']):
		for ii in parameters['imfs']:
			deltal = parameters['deltal_libs'][mi]
			model_wave_int,model_flux_int, age,metal = \
				get_model(parameters,mm,ii,deltal,data['vdisp'],data['wavelength'],data['r_instrument'],data['ebv_mw'])

			print "Matching data to models..."
			wave,data_flux,error_flux,model_flux_raw = \
				match_data_models(data['wavelength'],data['flux'],data['flags'],data['error'],model_wave_int,model_flux_int,\
					parameters['wave_limits'][0],parameters['wave_limits'][1])
			print "Normalising all spectra."
			model_flux,mass_factors = normalise_spec(data_flux,model_flux_raw)

			# Get filtered values IF dust is on!
			if parameters['hpf_mode']=='on':

				print "Determining attenuation curve through HPF fitting"
				best_ebv,attenuation_curve = determine_attenuation(wave,data_flux,error_flux,model_flux,parameters,age,metal)

				if parameters['plot_diagnostics']:
					print "Best ebv is "+str(best_ebv)
					plt.plot(attenuation_curve)
					plt.title("Attenuation curve")
					plt.show()


				# Apply curve to models and renormalise:
				print "Curve found! Applying to models..."
				model_flux_atten = np.zeros(np.shape(model_flux_raw))
				for m in range(len(model_flux_raw)):
					model_flux_atten[m] 	= attenuation_curve * model_flux_raw[m]
				model_flux,mass_factors = normalise_spec(data_flux,model_flux_atten)
				print "Fitting with attenuated models..."
				light_weights_int,chis_int,branch = fitter(wave,data_flux,error_flux,model_flux,parameters)
			
			elif parameters['hpf_mode'] == 'hpf_only':

				print "Using filtered values to determing SP properties only."
				smoothing_length = parameters['dust_smoothing_length']
				hpf_data    = hpf(data_flux)
				hpf_models  = np.zeros(np.shape(model_flux))
				for m in range(len(model_flux)):
					hpf_models[m] = hpf(model_flux[m])


				zero_dat = np.where( (np.isnan(hpf_data)) & (np.isinf(hpf_data)) )
				hpf_data[zero_dat] = 0.0
				for m in range(len(model_flux)):
					hpf_models[m,zero_dat] = 0.0
				hpf_error    = np.zeros(len(error_flux))
				hpf_error[:] = np.median(error_flux)/np.median(data_flux) * np.median(hpf_data)
				hpf_error[zero_dat] = np.max(hpf_error)*999999.9

				best_ebv = -99
				hpf_models,mass_factors = normalise_spec(hpf_data,hpf_models)
				light_weights_int,chis_int,branch = fitter(wave,hpf_data,hpf_error,hpf_models,parameters)

			elif parameters['hpf_mode'] == 'off':
				raise NotImplementedError("Not using a HPF and fitting using model curves not implemented yet")
				# use loop over dust curve, but this will take a while!
			
			

			print "Fitting complete! Calculating average properties and outputting."
			# Convert chis into probs
			# Degrees of freedom approximately = number of wavelength points
			dof = len(wave)
			probs_int 		= convert_chis_to_probs(chis_int,dof)
			# Remove zero-prob solutions
			nonzero_prob 	= np.where(probs_int>0.00001)
			probs 			= probs_int[nonzero_prob]
			light_weights 	= light_weights_int[nonzero_prob]
			chis 			= chis_int[nonzero_prob]

			# Get mass-weighted SSP contributions using saved M/L ratio. (raw and normalised)
			unnorm_mass,mass_weights = light_weights_to_mass(light_weights,mass_factors)
			

			# Calculate all average properties and errors
			averages = calculate_averages_pdf(probs,light_weights,mass_weights,unnorm_mass,\
												age,metal,parameters['pdf_sampling'],data['redshift'])
					
										

			unique_ages 				= np.unique(age)
			marginalised_age_weights 	= np.zeros(np.shape(unique_ages))

			marginalised_age_weights_int = np.sum(mass_weights.T,1)
			for ua in range(len(unique_ages)):
				marginalised_age_weights[ua] = np.sum(marginalised_age_weights_int[np.where(age==unique_ages[ua])])

			
			# sfr_int,sfr_error_int = star_formation_rate(np.log10(unique_ages)+9.0,marginalised_age_weights) 
			# sfr = sfr_int * 10**averages['stellar_mass'] / (10.0**7)
			# sfr_error = sfr_error_int * 10**averages['stellar_mass'] / (10.0**7)

			# print "Star formation rate is (in M / yr) "+str(sfr)+" plus.minus "+str(sfr_error)
			# Tracer()()

			best_fit_index = [np.argmin(chis)]
			best_fit = np.dot(light_weights[best_fit_index],model_flux)[0]
			


			if parameters['plot_fits']:
				plt.plot(wave,data_flux,'k')
				plt.plot(wave,best_fit,'r',linewidth=1.0)
				out_plot_string = 'plots/fit.eps'
				plt.savefig(out_plot_string,format='eps',transparent=True)
				plt.close()
			if parameters['plot_diagnostics']:
				plt.plot(wave,data_flux,'k')
				plt.plot(wave,best_fit,'r',linewidth=1.0)
				out_plot_string = 'plots/fit.eps'
				plt.show()
				plt.close()

			import plotting


			fits = np.dot(light_weights,model_flux)

			#Tracer()()
			#plotting.plot_fits(wave,data_flux,fits,probs)
			#plotting.plot_sfh_contours(age,metal,light_weights,probs,title="Light-weighted properties")
			#plotting.plot_sfh_contours(age,metal,mass_weights,probs,title="Mass-weighted properties")
			#Tracer()()
			
			# Calculate the weighted average of SSPs for the secondary outputs and contour plots.
			if parameters['observation_type']=='ifu':
				file1=parameters['output_dir_prefix']+parameters['file_in'].split('/')[-1]+'/'
				file2=parameters['output_dir_prefix']+parameters['file_in'].split('/')[-1]+'/'+mm+'/'
				file3=parameters['output_dir_prefix']+parameters['file_in'].split('/')[-1]+'/'+mm+'/'+ii+'/'
			else:
				file1=parameters['output_dir_prefix']
				file2=parameters['output_dir_prefix']+mm+'/'
				file3=parameters['output_dir_prefix']+mm+'/'+ii+'/'


			if not os.path.exists(file1):
				os.makedirs(file1)
			if not os.path.exists(file2):
				os.makedirs(file2)
			if not os.path.exists(file3):
				os.makedirs(file3)

			parameters['output_file'] 	= 	parameters['output_dir_prefix']+\
											parameters['file_in'].split('/')[-1]+'/'+mm+'/'+ii+'/'



			if parameters['observation_type']=='ifu':
				f = open(parameters['output_file']+'bin'+str(int(parameters['bin_number']))+'_single.txt', 'wb')
				f.write("# x, y, bin_number, Light_age / log(Gyrs) [value, +error, -error] light [Z/H] [value, +error,-error],"+\
							"mass age / log(Gyrs) [value, +error, -error],"+\
							"mass [Z/H][value, +error, -error], E(B-V), stellar mass [value, +error,-error]\n")
				
				f.write(str(data['xpos'])+'\t'+str(data['ypos'])+'\t'+str(parameters['bin_number'])+'\t'+\
						str(averages['light_age'])+'\t'+str(averages['light_age_1_sig_plus'])+'\t'+str(averages['light_age_1_sig_minus'])+'\t'+\
						str(averages['light_metal'])+'\t'+str(averages['light_metal_1_sig_plus'])+'\t'+str(averages['light_metal_1_sig_minus'])+'\t'+\
						str(averages['mass_age'])+'\t'+str(averages['mass_age_1_sig_plus'])+'\t'+str(averages['mass_age_1_sig_minus'])+'\t'+\
						str(averages['mass_metal'])+'\t'+str(averages['mass_metal_1_sig_plus'])+'\t'+str(averages['mass_metal_1_sig_minus'])+'\t'+\
						str(best_ebv)+'\t'+\
						str(averages['stellar_mass'])+'\t'+str(averages['stellar_mass_1_sig_plus'])+'\t'+str(averages['stellar_mass_1_sig_minus'])+'\n')

				f.close()
				print "Combining ascii fits files..."
				files_present = os.listdir(parameters['output_file'])
				string=str(parameters['file_in']).replace('./data/manga/','')
				name=string.replace('-LOGCUBE_BIN-RADIAL-015.fits','')

				if np.size(files_present)>0:
					
					combine_files = open(parameters['output_file']+'/'+name+'-combined.txt', 'wb')
					combine_files.write("# x, y, bin_number, Light_age / log(Gyrs) [value, +error, -error] light [Z/H] [value, +error,-error],"+\
							"mass age / log(Gyrs) [value, +error, -error],"+\
							"mass [Z/H][value, +error, -error], E(B-V), stellar mass [value, +error,-error]\n")
					for o in files_present:
						try:
							a = o.split('_')[-1]
						except IndexError:
							continue
						if o.split('_')[-1] == 'single.txt':
							fits 		= np.loadtxt(parameters['output_file']+o, skiprows=1, unpack=True)

							combine_files.write(str(fits[0])+'\t'+str(fits[1])+'\t'+str(fits[2])+\
												'\t'+str(fits[3])+'\t'+str(fits[4])+'\t'+str(fits[5])+\
												'\t'+str(fits[6])+'\t'+str(fits[7])+'\t'+str(fits[8])+\
												'\t'+str(fits[9])+'\t'+str(fits[10])+'\t'+str(fits[11])+\
												'\t'+str(fits[12])+'\t'+str(fits[13])+'\t'+str(fits[14])+\
												'\t'+str(fits[15])+'\t'+str(fits[16])+'\t'+str(fits[17])+'\t'+str(fits[18])+'\n')


					combine_files.close()
			else:

				f = open(parameters['output_dir_prefix']+mm+'/'+ii+'/'+parameters['file_in'].split('/')[-1]+'.txt', 'wb')
				f.write("# Light_age / log(Gyrs) [value, +error, -error] light [Z/H] [value, +error,-error],"+\
							"mass age / log(Gyrs) [value, +error, -error],"+\
							"mass [Z/H][value, +error, -errpr], E(B-V), stellar mass [value, +error,-error]\n")
				
				f.write(str(averages['light_age'])+'\t'+str(averages['light_age_1_sig_plus'])+'\t'+str(averages['light_age_1_sig_minus'])+'\t'+\
						str(averages['light_metal'])+'\t'+str(averages['light_metal_1_sig_plus'])+'\t'+str(averages['light_metal_1_sig_minus'])+'\t'+\
						str(averages['mass_age'])+'\t'+str(averages['mass_age_1_sig_plus'])+'\t'+str(averages['mass_age_1_sig_minus'])+'\t'+\
						str(averages['mass_metal'])+'\t'+str(averages['mass_metal_1_sig_plus'])+'\t'+str(averages['mass_metal_1_sig_minus'])+'\t'+\
						str(best_ebv)+'\t'+\
						str(averages['stellar_mass'])+'\t'+str(averages['stellar_mass_1_sig_plus'])+'\t'+str(averages['stellar_mass_1_sig_minus'])+'\n')

				f.close()
			print "Wrote ASCII output to "+parameters['output_file']

			



			
			#fitter(wavelength_in,data_in,error_in,models_in)

			# if options['plot_fits']:
			# 	#from plotting import *
			# 	plot_simple_fit()
			# 	plot_contours()

			# output_results(output_dir)
