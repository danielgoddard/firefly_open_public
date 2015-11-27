#match_data_models.py
import numpy as np
import copy
from IPython.core.debugger import Tracer

def match_data_models(data_wave_int,data_flux_int,data_flags,error_flux_int,model_wave_int,model_flux_int,\
						min_wave_in,max_wave_in):

	"""
	Take data wavelength, flux and model wavelength, fluxes
	and interpolate onto the bisection of lowest sampled array.
	Insert the interpolated points into high-res array, 
	then integrate between them to get a matched resolution array.

	"""  

	def bisect_array(array):
		bisected_array = np.zeros(len(array)-1)
		for ai in range(len(bisected_array)):
			bisected_array[ai] = (array[ai]+array[ai+1])/2.0
		return bisected_array

	num_models = len(model_flux_int)
	# Find most resolved wavelength array.
	min_wave = np.max([np.min(data_wave_int[np.where(data_flags==1)]),np.min(model_wave_int),min_wave_in])
	max_wave = np.min([np.max(data_wave_int[np.where(data_flags==1)]),np.max(model_wave_int),max_wave_in])

	print np.min(data_wave_int[np.where(data_flags==1)])
	print np.min(model_wave_int)
	print min_wave_in
	print "--"
	print np.max(data_wave_int[np.where(data_flags==1)])
	print np.max(model_wave_int)
	print max_wave_in
	print "----"

	# Number of points within each range
	loc_model 	= np.array(( model_wave_int <= max_wave) & (model_wave_int >= min_wave))
	if np.sum(loc_model)==0:
		raise ValueError("The wavelength range input is below or above model wavelength coverage!")
	model_wave 	= model_wave_int[loc_model]
	num_mod  	= np.sum(loc_model)
	model_flux 	= np.zeros((num_models,num_mod))

	#Tracer()()
	for m in range(num_models):
		model_flux[m] = model_flux_int[m][loc_model]
	
	loc_data 	= np.array(( data_wave_int <= max_wave) & (data_wave_int >= min_wave)) 
	if np.sum(loc_data)==0:
		raise ValueError("The wavelength range input is below or above data wavelength coverage!")
	num_dat  	= np.sum(loc_data)
	data_wave 	= data_wave_int[loc_data]
	data_flux 	= data_flux_int[loc_data]
	error_flux 	= error_flux_int[loc_data]


	# The one with the highest density is the high-res spectrum.

	if num_mod >= num_dat:
		print "More model points than data points! Downgrading models..."


		# Model = high-res spectrum. Insert bisected data into model for boundaries.
		bisect_data = bisect_array(data_wave)+np.min(data_wave)*0.0000000001
		# Interpolate the bisected model at these points.
		

		
		matched_model = np.zeros((num_models,len(bisect_data)-1))


		for m in range(num_models):
			model_flux_bounds 	= np.interp(bisect_data, model_wave, model_flux[m])
			combined_wave_int 	= np.concatenate((model_wave,bisect_data))
			combined_flux_int 	= np.concatenate((model_flux[m],model_flux_bounds))
			sort_indices 		= np.argsort(combined_wave_int)

			combined_wave 		= np.sort(combined_wave_int)
			boundary_indices 	= np.searchsorted(combined_wave,bisect_data)
			combined_flux 		= combined_flux_int[sort_indices]

			len_combo = len(combined_flux)

			for l in range(len(boundary_indices)-1):
				if boundary_indices[l+1] >= len_combo:
					matched_model[m][l] = matched_model[m][l-1]
				else:
					matched_model[m][l] = np.trapz(combined_flux[boundary_indices[l]:boundary_indices[l+1]+1],\
											x=combined_wave[boundary_indices[l]:boundary_indices[l+1]+1]) / \
											(combined_wave[boundary_indices[l+1]]-combined_wave[boundary_indices[l]])

		matched_wave = data_wave[1:-1]
		matched_data = data_flux[1:-1]
		matched_error = error_flux[1:-1]

		
	else:
		print "More data points than model points! Downgrading data..."
		# Data = high-res spectrum. Insert bisected model into data for boundaries.
		bisect_model = bisect_array(model_wave)+np.min(model_wave)*0.0000000001
		# Interpolate the bisected model at these points.
		boundaries 	= np.searchsorted(data_wave,bisect_model)

		data_flux_bounds 	= np.interp(bisect_model, data_wave, data_flux)
		error_flux_bounds 	= np.interp(bisect_model, data_wave, error_flux)
		combined_wave_int 	= np.concatenate((data_wave,bisect_model))
		combined_flux_int 	= np.concatenate((data_flux,data_flux_bounds))
		combined_error_int 	= np.concatenate((error_flux,error_flux_bounds))
		sort_indices 		= np.argsort(combined_wave_int)
	
		combined_wave 		= np.sort(combined_wave_int)
		boundary_indices 	= np.searchsorted(combined_wave,bisect_model)
		combined_flux 		= combined_flux_int[sort_indices]
		combined_error 		= combined_error_int[sort_indices]

		matched_data,matched_error= np.zeros(len(boundary_indices)-1),np.zeros(len(boundary_indices)-1)

		len_combo = len(combined_flux)
		for l in range(len(boundary_indices)-1):
			if boundary_indices[l+1] >= len_combo:
				matched_data[l] 	= matched_data[l-1]
				matched_error[l] 	= matched_error[l-1]

			else:

				matched_data[l] 	= np.trapz(combined_flux[boundary_indices[l]:boundary_indices[l+1]+1], \
										x=combined_wave[boundary_indices[l]:boundary_indices[l+1]+1])/ \
											(combined_wave[boundary_indices[l+1]]-combined_wave[boundary_indices[l]])

				matched_error[l] 	= np.trapz(combined_error[boundary_indices[l]:boundary_indices[l+1]+1], \
											x=combined_wave[boundary_indices[l]:boundary_indices[l+1]+1])/ \
											(combined_wave[boundary_indices[l+1]]-combined_wave[boundary_indices[l]])

		matched_wave 		= model_wave[1:-1]
		matched_model 		= np.zeros((num_models,len(matched_wave)))
		for m in range(num_models):
			matched_model[m][:] = model_flux[m][1:-1]




	return matched_wave,matched_data,matched_error,matched_model





