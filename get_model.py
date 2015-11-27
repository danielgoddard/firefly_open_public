#get_model.py
import numpy as np
import pandas as pd
import astropy.io.fits as pyfits
import os
import glob
import copy
from idlutils import *

from IPython.core.debugger import Tracer
from downgrader import *
import matplotlib.pyplot as plt
import time
from dust import *

def get_model(options,model_used,imf_used,deltal,vdisp,wave_instrument,r_instrument,ebv_mw):

	"""
	Retrieves all relevant model files, in their downgraded format.
	If they aren't downgraded to the correct resolution / velocity dispersion,
	takes the base models in their native form and converts to downgraded files.

    If downgrading, needs the resolution profile of the instrument (wavelength,R).
	"""

	if options['models'] == 'm11':

		# Find the nearest 5km/s velocity dispersion.
		vdisp_round = int(round(vdisp/5.0)*5.0)
		#age_used 	= options['age_used']
		#metal_used 	= options['metal_used']

		# Define model fluxes once one run is done with this switch:
		first_file  = True 
		      
		
		model_files = []

		# FIND NUMBER OF METALLICITIES, AGES (metallicities in separate files)
		if model_used == 'MILES_UVextended' or model_used == 'MILES_revisedIRslope':
			model_path 		= 'models/SSP_M11_MILES/ssp_M11_'+model_used+'.'+imf_used
		else:
			model_path 		= 'models/SSP_M11_'+model_used+'/ssp_M11_'+model_used+'.'+imf_used
		# List all metallicities that are here:
		all_metal_files = glob.glob(model_path+'*')
		metal_files 	= []
		metal 			= []
		for z in range(len(all_metal_files)):
			zchar = all_metal_files[z][len(model_path):]
			if zchar == 'z001':
				znum = -0.3
			elif zchar == 'z002':
				znum = 0.0
			elif zchar == 'z004':
				znum = 0.3
			elif zchar == 'z0001.bhb':
				znum = -1.301
			elif zchar == 'z0001.rhb':
				znum = -1.302
			elif zchar == 'z10m4.bhb':
				znum = -2.301
			elif zchar == 'z10m4.rhb':
				znum = -2.302
			elif zchar == 'z10m4':
				znum = -2.300
			else:
				raise NameError('Unrecognised metallicity! Check model file names.')

			if znum>options['Z_limits'][0] and znum<options['Z_limits'][1]:
				metal_files.append(all_metal_files[z])
				metal.append(znum)



		model_flux, age_model, metal_model = [],[],[]

		print "Will redden models by the F99 curve with E(B-V) = "+str(ebv_mw)
		for zi,z in enumerate(metal_files):
			print "Retrieving and downgrading models for "+z
			model_table = pd.read_table(z,converters={'Age':np.float64},header=None,usecols=[0,2,3],\
						names=['Age','wavelength_model','flux_model'],delim_whitespace=True)
			# GET WAVELENGTH AND FLUX OF MODELS

			# Get list of ages:
			age_data = np.unique(model_table['Age'].values.ravel())

			for a in age_data:
				logyrs_a = np.log10(a)+9.0
				if logyrs_a < options['age_limits'][0] or logyrs_a > options['age_limits'][1]:
					continue
				spectrum 		= model_table.loc[model_table.Age==a,['wavelength_model','flux_model']].values
				#wavelength,flux = spectrum[:,0],spectrum[:,1]
				wavelength_int,flux = spectrum[:,0],spectrum[:,1]

				if options['data_wave_medium'] == 'vacuum':
					wavelength = airtovac(wavelength_int)
				else:
					wavelength = wavelength_int

				
				mf = downgrade(wavelength,flux,deltal,vdisp_round,wave_instrument,r_instrument)
				
				#string =z.replace('models/SSP_M11_MILES/','')
				#f, axarr = plt.subplots(2, sharex=True)
				#axarr[0].plot(wavelength, flux)
				#axarr[1].plot(wavelength, mf)
				#plt.savefig('/Users/Daniel/firefly_original/downgraded_plots_original/'+string+'.png')
				#plt.close()

			
				
				# Supply negative ebv_mw to redden the models
				if ebv_mw != 0:

					attenuations = unred(wavelength,ebv=0.0-ebv_mw)
					# if options['plot_diagnostics']:
					# 	plt.plot(mf)
					# 	plt.plot(mf*attenuations)
					# 	plt.show()  
					model_flux.append(mf*attenuations)
				else:
					model_flux.append(mf)

				age_model.append(a)
				metal_model.append(metal[zi])
				first_model = False

			#print str(mf)+" of "+str(num_models)

		# if first_file:
		# 	raise NameError('Cannot find model files (base or downgraded). Stopping.')

		print "Retrieved all models!"
		if np.size(metal_files) == 0:
			raise ValueError("No models with this set of parameters available!")
		return wavelength,model_flux,age_model,metal_model

	elif options['models'] == 'm09':

		# Find the nearest 5km/s velocity dispersion.
		vdisp_round = int(round(vdisp/5.0)*5.0)
		#age_used 	= options['age_used']
		#metal_used 	= options['metal_used']

		# Define model fluxes once one run is done with this switch:
		first_file  = True 
		      
		
		model_files = []

		# FIND NUMBER OF METALLICITIES, AGES (metallicities in separate files)
		if options['downgrade_models']:
			model_path 		= 'models/UVmodels_Marastonetal08b/'
		else:
			model_path 		= 'models/UVmodels_Marastonetal08b_downgraded/'
		# List all metallicities that are here:
		all_metal_files = glob.glob(model_path+'*')
		metal_files 	= []
		metal 			= []
		for z in range(len(all_metal_files)):
			zchar = all_metal_files[z].split('.')[1][2:]
			if zchar == 'z001':
				znum = -0.3
			elif zchar == 'z002':
				znum = 0.0
			elif zchar == 'z004':
				znum = 0.3
			elif zchar == 'z0001':
				znum = -1.300
			else:
				raise NameError('Unrecognised metallicity! Check model file names.')

			if znum>options['Z_limits'][0] and znum<options['Z_limits'][1]:
				metal_files.append(all_metal_files[z])
				metal.append(znum)



		model_flux, age_model, metal_model = [],[],[]

		print "Will redden models by the F99 curve with E(B-V) = "+str(ebv_mw)
		for zi,z in enumerate(metal_files):
			print "Retrieving and downgrading models for "+z
			model_table = pd.read_table(z,converters={'Age':np.float64},header=None,usecols=[0,1,2],\
						names=['Age','wavelength_model','flux_model'],delim_whitespace=True)
			# GET WAVELENGTH AND FLUX OF MODELS

			# Get list of ages:
			age_data = np.unique(model_table['Age'].values.ravel())

			for a in age_data:
				logyrs_a = np.log10(a)+9.0
				if logyrs_a < options['age_limits'][0] or logyrs_a > options['age_limits'][1]:
					continue
				spectrum 		= model_table.loc[model_table.Age==a,['wavelength_model','flux_model']].values
				wavelength_int,flux = spectrum[:,0],spectrum[:,1]

				if options['data_wave_medium'] == 'vacuum':
					wavelength = airtovac(wavelength_int)
				else:
					wavelength = wavelength_int

				if options['downgrade_models']:
					mf = downgrade(wavelength,flux,deltal,vdisp_round,wave_instrument,r_instrument)
				else:
					mf = copy.copy(flux)			

				# Supply negative ebv_mw to redden the models
				if ebv_mw != 0:

					attenuations = unred(wavelength,ebv=0.0-ebv_mw)
					# if options['plot_diagnostics']:
					# 	plt.plot(mf)
					# 	plt.plot(mf*attenuations)
					# 	plt.show()  
					model_flux.append(mf*attenuations)
				else:
					model_flux.append(mf)

				age_model.append(a)
				metal_model.append(metal[zi])
				first_model = False

			#print str(mf)+" of "+str(num_models)

		# if first_file:
		# 	raise NameError('Cannot find model files (base or downgraded). Stopping.')

		print "Retrieved all models!"
		return wavelength,model_flux,age_model,metal_model

	elif options['models']=='bc03':

		# Find the nearest 5km/s velocity dispersion.
		vdisp_round = int(round(vdisp/5.0)*5.0)

		# Define model fluxes once one run is done with this switch:
		first_file  = True 
		      
		
		model_files = []

		# FIND NUMBER OF METALLICITIES, AGES (metallicities in separate files)
		model_path 		= 'models/spectra_bc03/'
		# List all metallicities that are here:
		all_metal_files = glob.glob(model_path+'*')
		metal_files 	= []
		metal_model 	= []
		age_model 		= []
		model_flux 		= []

		for z in range(len(all_metal_files)):
			if all_metal_files[z][0]=='.':
				continue

			junk1,junk2,modelchar,imfchar,zchar,agechar = all_metal_files[z].split('/')[-1].split('_')
			if zchar == 'Z002':
				znum = 0.0
			elif zchar == 'Z0004':
				znum = np.log10(4.0/20.0)
			elif zchar == 'Z005':
				znum = np.log10(2.5)
			elif zchar == 'Z0008':
				znum = np.log10(8.0/20.0)
			elif zchar == 'z10m4.rhb':
				znum = -2.302
			else:
				raise NameError('Unrecognised metallicity! Check model file names.')

			age_num = np.log10(float(agechar[3:-4]))+9.0

			if znum>options['Z_limits'][0] and znum<options['Z_limits'][1] and \
				age_num>options['age_limits'][0] and age_num<options['age_limits'][1]:
				metal_files.append(all_metal_files[z])
				metal_model.append(znum)
				age_model.append(age_num)




		print "Will redden models by the F99 curve with E(B-V) = "+str(ebv_mw)
		for zi,z in enumerate(metal_files):
			print "Retrieving and downgrading models for "+z
			model_table = pd.read_table(z,header=None,usecols=[0,1],skiprows=6,\
						names=['wavelength_model','flux_model'],delim_whitespace=True)
			# GET WAVELENGTH AND FLUX OF MODELS

			spectrum 		= model_table.values
			wavelength_int,flux = spectrum[400:-1000,0],spectrum[400:-1000,1]/(3.826*10.0**33)

			if options['data_wave_medium'] == 'vacuum':
				wavelength = airtovac(wavelength_int)
			else:
				wavelength = wavelength_int
			mf = downgrade(wavelength,flux,deltal,vdisp_round,wave_instrument,r_instrument)
			

			# Supply negative ebv_mw to redden the models
			if ebv_mw != 0:

				attenuations = unred(wavelength,ebv=0.0-ebv_mw)
				# if options['plot_diagnostics']:
				# 	plt.plot(mf)
				# 	plt.plot(mf*attenuations)
				# 	plt.show()  
				model_flux.append(mf*attenuations)
			else:
				model_flux.append(mf)

			
			#print str(mf)+" of "+str(num_models)

		# if first_file:
		# 	raise NameError('Cannot find model files (base or downgraded). Stopping.')

		print "Retrieved all models!"
		return wavelength,model_flux,age_model,metal_model


	else:
		"Other models are not implemented yet."
