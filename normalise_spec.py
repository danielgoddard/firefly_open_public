#normalise_spec.py
import numpy as np

def normalise_spec(data_flux,model_flux):

	"""
	Normalises all spectra (already clipped) to their median values.
	Saves the factors for later use.
	"""

	data_norm 				= np.median(data_flux)
	num_mods 				= len(model_flux)
	model_norm,mass_factor 	= np.zeros(num_mods),np.zeros(num_mods)
	normed_model_flux 		= np.zeros((num_mods,len(model_flux[0])))

	for m in range(len(model_flux)):
		model_norm[m] 			= np.median(model_flux[m])
		mass_factor[m] 			= data_norm/model_norm[m]
		normed_model_flux[m] 	= model_flux[m] / model_norm[m] * data_norm

	return normed_model_flux,mass_factor