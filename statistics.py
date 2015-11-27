import numpy as np
from scipy.stats import chi2 
import matplotlib.pyplot as plt
from astropy.cosmology import Planck13
from IPython.core.debugger import Tracer
import math
import matplotlib.pyplot as plt
import statistics
import plotting

def convert_chis_to_probs(chis,dof):

	chis = chis / np.min(chis) * dof
	prob =  1.0 - chi2.cdf(chis,dof)
	prob = prob / np.max(prob)

	return prob

def light_weights_to_mass(light_weights,mass_factors):
	"""
	Uses the data/model mass-to-light ratio to convert
	SSP contribution (weights) by light into 
	SSP contributions by mass.
	"""
	
	mass_weights 	= np.zeros(np.shape(light_weights))
	unnorm_mass 	= np.zeros(np.shape(light_weights))

	for w in range(len(light_weights)):
		# Check this is the correct way around:
		unnorm_mass[w]	= light_weights[w] * mass_factors
		mass_weights[w] = unnorm_mass[w] / np.sum(unnorm_mass[w])

	return unnorm_mass,mass_weights


def calculate_averages_pdf(probs,light_weights,mass_weights,unnorm_mass,age,metal,sampling,redshift):

	"""
	Calculates light- and mass-averaged age and metallicities.
	Also outputs stellar mass and mass-to-light ratios.
	And errors on all of these properties.

	It works by taking the complete set of probs-properties and
	maximising over the parameter range (such that solutions with
	equivalent values but poorer probabilities are excluded). Then,
	we calculate the median and 1/2 sigma confidence intervals from 
	the derived 'max-pdf'.

	NB: Solutions with identical SSP component contributions 
	are re-scaled such that the sum of probabilities with that
	component = the maximum of the probabilities with that component.
	i.e. prob_age_ssp1 = max(all prob_age_ssp1) / sum(all prob_age_ssp1) 
	This is so multiple similar solutions do not count multiple times.

	Outputs a dictionary of:
	- light_[property], light_[property]_[1/2/3]_sigerror
	- mass_[property], mass_[property]_[1/2/3]_sigerror
	- stellar_mass, stellar_mass_[1/2/3]_sigerror
	- mass_to_light, mass_to_light_[1/2/3]_sigerror
	- maxpdf_[property]
	- maxpdf_stellar_mass
	where [property] = [age] or [metal]
	"""

	# Sampling number of max_pdf (100:recommended) from options
	

	def averages_and_errors(probs,prop,sampling):
		def max_pdf(probs,prop,sampling):

			lower_limit 	= np.min(prop)
			upper_limit 	= np.max(prop)
			if upper_limit==lower_limit:
				return np.asarray(prop),np.ones(len(probs))/np.size(probs)
			property_pdf_int= np.arange(lower_limit,upper_limit*1.001,(upper_limit-lower_limit)/sampling)+(upper_limit-lower_limit)*0.00000001
			
			prob_pdf 		= np.zeros(len(property_pdf_int))

			def bisect_array(array):
				bisected_array = np.zeros(len(array)-1)
				for ai in range(len(bisected_array)):
					bisected_array[ai] = (array[ai]+array[ai+1])/2.0
				return bisected_array


			for p in range(len(property_pdf_int)-1):
				match_prop = np.where( (prop <= property_pdf_int[p+1]) & (prop > property_pdf_int[p]) )
				if np.size(match_prop) == 0:
					continue
				else:
					prob_pdf[p] = np.max( probs[match_prop] )

			property_pdf = bisect_array(property_pdf_int)
			return property_pdf,prob_pdf[:-1]/np.sum(prob_pdf)
		
		print' we are now in the averages section.'
		
		best_fit 		= np.argmax(probs)
		one_sig_fits	= np.where(probs>(1-0.6827))
		weights=light_weights
		fits_arr = [best_fit]
		for i,fits in enumerate(fits_arr):

			if i == 0:
				ii = (0,0)
			
			inc_weights = statistics.renormalise_weights(weights[fits_arr],probs[fits_arr])
			inc_weights = inc_weights / np.sum(inc_weights)
		
			prev_weight = np.zeros(len(weights[0]))
			tot_weight 	= np.sum(inc_weights,0)
				
		
		print sampling
		print '-'
		
		xdf,y = max_pdf(probs,prop,sampling)
		cdf = np.zeros(np.shape(y))
		cdf_probspace = np.zeros(np.shape(y))



		for m in range(len(y)):
			cdf[m] = np.sum(y[:m])
		cdf = cdf / np.max(cdf)

		#fig,ax=plt.subplots()
		#ax.tick_params(labelsize=20)
		#ax.plot(prop+9.0,probs/np.max(probs),'k.',label='Individual fits')
		#ax.plot(xdf+9.0,y/np.max(y),'r--',linewidth=2.0,label='Likelihood distribution')

		#xlim_in = (np.min(xdf[np.where(y/np.max(y)>0.01)])+9.0,np.max(xdf[np.where(y/np.max(y)>0.01)])+9.0)
		#plt.plot(xdf,cdf,'b-')
		#ax.set_ylabel("Likelihood",fontsize=22)
		#ax.set_xlabel("Age / log(yrs)",fontsize=22)
		#ax.legend(loc='upper left', shadow=True,prop={'size':20})
		#out_plot_string = './plots/likelihood_'+str(np.mean(prop))+'.eps'
		#plt.savefig(out_plot_string,format='eps',transparent=True,bbox_inches='tight')
		#plt.close()

		area_probspace = y*(xdf[1]-xdf[0])
		area_probspace = area_probspace/np.sum(area_probspace)


		indx_probspace = np.argsort(area_probspace)[::-1]
		desc_probspace = np.sort(area_probspace)[::-1]

		cdf_probspace = np.zeros(np.shape(desc_probspace))
		for m in range(len(desc_probspace)):
			cdf_probspace[m] = np.sum(desc_probspace[:m])

		def find_closest(A, target):
			#A must be sorted
			idx = A.searchsorted(target)
			idx = np.clip(idx, 1, len(A)-1)
			left = A[idx-1]
			right = A[idx]
			idx -= target - left < right - target
			return idx

		# Median, + / - 1 sig, + / - 2 sig, + / - 3 sig
		#av_sigs = [0.5,0.6827/2.0+0.5,0.5-0.6827/2.0,0.95/2.0+0.5,0.5-0.95/2.0,0.997/2.0+0.5,0.5-0.997/2.0]
		av_sigs = [0.6827,0.9545,0.9973]

		
		# Sorts results by likelihood and calculates confidence intervals on sorted space
		index_close = find_closest(cdf_probspace, av_sigs)
		
		best_fit 					= xdf[indx_probspace[0]]
		upper_onesig,lower_onesig 	= np.max(xdf[indx_probspace[:index_close[0]]]),np.min(xdf[indx_probspace[:index_close[0]]])
		upper_twosig,lower_twosig 	= np.max(xdf[indx_probspace[:index_close[1]]]),np.min(xdf[indx_probspace[:index_close[1]]])
		upper_thrsig,lower_thrsig 	= np.max(xdf[indx_probspace[:index_close[2]]]),np.min(xdf[indx_probspace[:index_close[2]]])

		# Takes whole pdf and computes median, confidence intervals
		#index_close = find_closest(cdf, av_sigs)


		if np.size(xdf) == 0:
			raise Exception('No solutions found??? FIREFLY error (see statistics.py)')
		
		#return xdf[index_close]
		return [best_fit,upper_onesig,lower_onesig,upper_twosig,lower_twosig,upper_thrsig,lower_thrsig]


	log_age = np.log10(age)
	log_age[np.isnan(log_age)|np.isinf(log_age)] = 0.0

	av = {}

	av['light_age'],av['light_age_1_sig_plus'],av['light_age_1_sig_minus'],\
		 	  av['light_age_2_sig_plus'],av['light_age_2_sig_minus'],\
		 	  av['light_age_3_sig_plus'],av['light_age_3_sig_minus'] \
				= averages_and_errors(probs,np.dot(light_weights,log_age),sampling)

	av['light_metal'],av['light_metal_1_sig_plus'],av['light_metal_1_sig_minus'],\
		 	  av['light_metal_2_sig_plus'],av['light_metal_2_sig_minus'],\
		 	  av['light_metal_3_sig_plus'],av['light_metal_3_sig_minus'] \
				= averages_and_errors(probs,np.dot(light_weights,metal),sampling)

	av['mass_age'],av['mass_age_1_sig_plus'],av['mass_age_1_sig_minus'],\
		 	 av['mass_age_2_sig_plus'],av['mass_age_2_sig_minus'],\
		 	 av['mass_age_3_sig_plus'],av['mass_age_3_sig_minus'] \
				= averages_and_errors(probs,np.dot(mass_weights,log_age),sampling)

	av['mass_metal'],av['mass_metal_1_sig_plus'],av['mass_metal_1_sig_minus'],\
		 	 av['mass_metal_2_sig_plus'],av['mass_metal_2_sig_minus'],\
		 	 av['mass_metal_3_sig_plus'],av['mass_metal_3_sig_minus'] \
				= averages_and_errors(probs,np.dot(mass_weights,metal),sampling)

	dist_lum 			= Planck13.luminosity_distance(redshift)*10.0**6 * 3.09 * 10.0**18
	conversion_factor 	= 10.0**(-17) * 4*math.pi*dist_lum.value**2.0

	tot_mass = np.log10(np.sum(unnorm_mass,1) * conversion_factor)
	av['stellar_mass'],av['stellar_mass_1_sig_plus'],av['stellar_mass_1_sig_minus'],\
		 	 av['stellar_mass_2_sig_plus'],av['stellar_mass_2_sig_minus'],\
		 	 av['stellar_mass_3_sig_plus'],av['stellar_mass_3_sig_minus'] \
				= averages_and_errors(probs,tot_mass,sampling)
	
	#plotting.files_of_light_weights(log_age,metal,light_weights,probs,"output_test")
	#plotting.files_of_mass_weights(log_age,metal,light_weights,probs,"output_test_mass")
	
	return av


def renormalise_weights(weights,probs):

	# Multiple solutions with same base can contribute to a fit multiply. Simply
	# re-scale these solutions so that their weights is the total of
	# all matching solution probabilities divided by their max probability..
	weights_base 	= weights.T
	probs_factor 	= np.zeros(len(weights_base))
	renorm_weights 	= np.zeros(np.shape(weights))

	for w in range(len(weights_base)): 
		match_index = np.where(weights_base[w]>0)
		if np.size(match_index)==0:
			probs_factor[w]=1
		else:
			probs_factor[w] = np.max(probs[match_index])/np.sum(probs[match_index])

	for f in range(len(probs)):
		renorm_weights[f] = weights[f]*np.prod(probs_factor[np.nonzero(weights[f])])

	return renorm_weights








