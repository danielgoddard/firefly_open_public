import numpy as np
import matplotlib.pyplot as plt
#import seaborn
import scipy
import scipy.ndimage
import matplotlib
from matplotlib.cm import get_cmap
from scipy.interpolate import griddata
import matplotlib.colors as colors
from matplotlib import gridspec
from matplotlib import rc
import statistics

"""
A selection of plotting tools used in FIREFLY's standard outputs.
"""


def plot_simple():
	# Plot a 1D spectrum.
	pass

def plot_simple_fit():
	# Plot a data, model, and error spectral fit.
	pass

def plot_fits(wave,data,fits,probs,title="",min_wave=0,max_wave=999999999):
	# Plot the spectral fits to the data with probability weighted fit 
	# shown with 1 sigma confidence intervals.
	# Input: model fits, data, error, probabilities.

	# 
	if min_wave!=0 or max_wave!=999999999:

		allowed_wave = np.where( ( wave>min_wave) & (wave <max_wave) )[0]
		wave = wave[allowed_wave]
		data = data[allowed_wave]

		fits_new = np.zeros((len(fits),len(allowed_wave)))
		for f in range(len(fits)):
			fits_new[f] = fits[f][allowed_wave]
		fits = fits_new

	fig, (ax,ax_residuals) = plt.subplots(2, sharex=True)
	gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1]) 
	ax = plt.subplot(gs[0])
	ax_residuals = plt.subplot(gs[1])
	ax.plot(wave,data,'k-',linewidth=0.5,label="Data")

	ax.tick_params(labelsize=16)
	ax.set_xticklabels("")
	ax_residuals.tick_params(labelsize=16)
	ax.set_ylabel(r"Flux / erg s$^{-1}$cm$^{-2}$$\AA^{-1}$",fontsize=20)

	ax_residuals.set_xlabel(r"Wavelength / $\AA^{-1}$",fontsize=20)
	ax_residuals.set_ylabel("Model / data",fontsize=18)

	best_fit 		= np.argmax(probs)
	one_sig_fits	= np.where(probs>(1-0.6827))
	two_sig_fits 	= np.where(probs>(1-0.9545))
	thr_sig_fits 	= np.where(probs>(1-0.997))

	ax.plot(wave,fits[best_fit],'r',linewidth=.5,label="Best fit")

	ax.set_xlim(np.min(wave),np.max(wave))
	ax_residuals.set_xlim(np.min(wave),np.max(wave))
	ax.set_ylim(np.min(data)*0.95,np.max(data)*1.05)
	upper_onesig = np.max(fits[one_sig_fits],0)
	lower_onesig = np.min(fits[one_sig_fits],0)

	upper_twosig = np.max(fits[two_sig_fits],0)
	lower_twosig = np.min(fits[two_sig_fits],0)

	upper_thrsig = np.max(fits[thr_sig_fits],0)
	lower_thrsig = np.min(fits[thr_sig_fits],0)

	ax.plot([0,0],[0,0],'-',color='green', alpha=0.5,label="68% level",linewidth=10)
	ax.plot([0,0],[0,0],'-',color='blue', alpha=0.5,label="95% level",linewidth=10)
	ax.fill_between(wave, lower_onesig, upper_onesig, facecolor='green', alpha=0.5,label="68% level")
	ax.fill_between(wave, lower_onesig, lower_twosig, facecolor='blue', alpha=0.5,label="95% level")
	ax.fill_between(wave, upper_onesig, upper_twosig, facecolor='blue', alpha=0.5)

	ax_residuals.plot([np.min(wave),np.max(wave)],[0,0],'k-',linewidth=0.5)
	ax_residuals.plot(wave,fits[best_fit]/data,'r-',linewidth=0.5)
	ax_residuals.fill_between(wave, lower_onesig/data, upper_onesig/data, facecolor='green', alpha=0.5)
	ax_residuals.fill_between(wave, lower_onesig/data, lower_twosig/data, facecolor='blue', alpha=0.5)
	ax_residuals.fill_between(wave, upper_onesig/data, upper_twosig/data, facecolor='blue', alpha=0.5)

	ax.legend(loc=[0.4,0.05],ncol=2,prop={'size':16},shadow=True)

	ax_residuals.set_ylim(0,2)

	# ax.fill_between(wave, lower_twosig, lower_thrsig, facecolor='beige', alpha=0.5)
	# ax.fill_between(wave, upper_twosig, upper_thrsig, facecolor='beige', alpha=0.5)

	fig.subplots_adjust(hspace=0)

	out_plot_string = 'plots/fits_'+title+'.eps'
	plt.savefig(out_plot_string,format='eps',transparent=True,bbox_inches='tight')
	plt.close()


def plot_sfh_contours(age,metal,weights,probs,title=""):
	# Plot FIREFLY'S likelihood-weighted SSP contribution contours
	# INPUT: ages, metallicities, weights, probs

	age = np.log10(age)+9.0
	metal = np.asarray(metal)
	def diff_array(array):
		bisected_array = np.zeros(len(array)-1)
		for ai in range(len(bisected_array)):
			bisected_array[ai] = (array[ai+1]-array[ai])
		return bisected_array

	min_width = np.min(np.absolute((diff_array(age))))

	best_fit 		= np.argmax(probs)
	one_sig_fits	= np.where(probs>(1-0.6827))
	two_sig_fits 	= np.where(probs>(1-0.9545))


	fits_arr = [[best_fit],one_sig_fits,two_sig_fits]
	label_arr = ["Best fit","68% level","95% level"]
	col_cum = ['black','green','cyan']
	sort_ages 	= np.sort(age)
	sort_age_index = np.argsort(age)
	sort_unique_ages = np.sort(np.unique(age))
	sort_metals = np.sort(metal)
	sort_unique_metals = np.sort(np.unique(metal))
	name_metals = [str(x) for x in sort_unique_metals]

	f, ax_arr = plt.subplots(2, 2, sharex=True, sharey=True)
	# ax_arr = ( (1,2),(3,4) )
	
	def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
		new_cmap = colors.LinearSegmentedColormap.from_list(
		'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
		cmap(np.linspace(minval, maxval, n)))
		return new_cmap

	cmap=truncate_colormap(get_cmap('hot_r'),minval=0.15,maxval=1.0,n=1000)

	cmap.set_under('w')
	bounds_interp_age = (9.0,10.19)
	bounds_interp_metal = (-2.5,0.90)

	ax_arr[0,1]= plt.subplot(222)
	ax_arr[0,1].set_xlim(bounds_interp_age[0],bounds_interp_age[1])
	#ax_arr[0,1].set_xlabel("Age / log(yr)",fontsize=20)

	for i,fits in enumerate(fits_arr):

		#ii = (int(i/2),i%2)

		if i == 0:
			ii = (0,0)
		elif i == 1:
			ii = (1,0)
		elif i ==2:
			ii = (1,1)

		ax_arr[ii].tick_params(labelsize=14)
		ax_arr[ii].plot([0,0],[0,0],'-',color=col_cum[i],markersize=1,label=label_arr[i],linewidth=3.0)
		#inc_weights = weights[fits]
		import statistics
		inc_weights = statistics.renormalise_weights(weights[fits],probs[fits])
		inc_weights = inc_weights / np.sum(inc_weights)



		prev_weight = np.zeros(len(weights[0]))
		tot_weight 	= np.sum(inc_weights,0)

		index_allow = np.where( (age>0) & (metal>-90) )
		npts 	= 300

		

		extent_array = [bounds_interp_age[0],bounds_interp_age[1],bounds_interp_metal[0],bounds_interp_metal[1]]

		xage 	= np.linspace(bounds_interp_age[0],bounds_interp_age[1],npts)
		ymetal 	= np.linspace(bounds_interp_metal[0],bounds_interp_metal[1],npts)
		ax_arr[ii].set_xlim(bounds_interp_age[0],bounds_interp_age[1])
		ax_arr[ii].set_ylim(bounds_interp_metal[0],bounds_interp_metal[1])
		# grid the data.
		zi = griddata((age, metal), tot_weight, (xage[None,:], ymetal[:,None]), method='linear')
		zi[np.isnan(zi)|np.isinf(zi)]=0.00001
		zs =scipy.ndimage.gaussian_filter(zi,sigma=5,order=0)
		# im = ax_arr[i].imshow(zi.T[::-1], origin='lower',interpolation='sinc', cmap=cmap,extent=extent_array,\
		# 						vmin=0,vmax=1,aspect='auto')
		print np.max(zi)
		CS = ax_arr[ii].contour(xage,ymetal,zs/np.max(zs),[0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],\
								linewidths=0.5,colors='k',vmin=0.01)
		CS = ax_arr[ii].contourf(xage,ymetal,zs/np.max(zs),100,cmap=cmap,vmin=0.01,vmax=1.0)
		ax_arr[ii].legend(loc="lower left",prop={'size':18},shadow=False)


		# Cumulative plots
		age_pdf = np.sum(zi,0)
		cum_pdf = np.zeros(np.shape(age_pdf))
		for ap in range(len(age_pdf)):
			cum_pdf[ap] = np.sum(age_pdf[:ap])


		ax_arr[0,1].plot(xage,cum_pdf/np.max(cum_pdf),'-',linewidth=3.0,color=col_cum[i]) 
		

	ax_arr[0,0].set_title(title,fontsize=20,y=1.05,x=1.0)

	ax_arr[0,1].plot([0,0],[0,0],'w-',markersize=0.1,label='CDF',linewidth=3.0)
	ax_arr[0,1].legend(loc="upper left",prop={'size':18},shadow=True,handlelength=-0.0,handletextpad=0.0)
	ax_arr[0,1].set_xticklabels([""])
	ax_arr[0,1].tick_params(labelsize=14)
	ax_arr[0,1].set_ylabel("Normalised Weight (< Age)",fontsize=16,rotation=270,labelpad=30)
	ax_arr[0,1].yaxis.set_label_position("right")
	ax_arr[0,1].yaxis.tick_right()

	ax_arr[1,0].set_ylabel("[Z/H]",fontsize=20,y=1.0)
	ax_arr[1,0].set_xlabel("Age / log(yr)",fontsize=20,x=1.0)
	plt.tight_layout()
	f.subplots_adjust(hspace=0,wspace=0)
	f.subplots_adjust(bottom=0.2)
	cbar_ax = f.add_axes([0.10,0.03,0.8,0.03])#[0.05, 0.15, 0.03, 0.8])
	cb = f.colorbar(CS, cax=cbar_ax,cmap=cmap,orientation='horizontal')
	cb.cmap.set_under(color='white', alpha=None)
	cb.set_ticks([0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
	cb.set_ticklabels(['0.01','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'])
	cbar_ax.tick_params(length=13, width=1, which='major')
	cb.set_label("Relative SSP weight",fontsize=20,labelpad=5)
	cb.ax.tick_params(labelsize=18)
	out_plot_string = 'plots/contour_'+title+'.eps'
	plt.savefig(out_plot_string,format='eps',transparent=True,bbox_inches='tight')
	plt.close()
		# Make a CDF:


def plot_sfh(age,metal,weights,probs):
	# Plot histograms of the SFH of the best fit, 1 sig, (2 sig?) solutions.
	# (Maybe with CDF)
	# INPUT: ages, metallicities, weights, probs

	age = np.log10(age)+9.0
	def diff_array(array):
		bisected_array = np.zeros(len(array)-1)
		for ai in range(len(bisected_array)):
			bisected_array[ai] = (array[ai+1]-array[ai])
		return bisected_array

	min_width = np.min(np.absolute((diff_array(age))))

	best_fit 		= np.argmax(probs)
	one_sig_fits	= np.where(probs>(1-0.6827))
	two_sig_fits 	= np.where(probs>(1-0.9545))

	fits_arr = [[best_fit],one_sig_fits,two_sig_fits]
	print fits_arr
	sort_ages 	= np.sort(age)
	sort_age_index = np.argsort(age)
	sort_unique_ages = np.sort(np.unique(age))
	sort_metals = np.sort(metal)
	sort_unique_metals = np.sort(np.unique(metal))
	color_metals = ['violet','blue','MediumAquaMarine','green','yellow','orange','darkred']
	name_metals = [str(x) for x in sort_unique_metals]

	f, ax_arr = plt.subplots(3, sharex=True, sharey=False)
	ax_cdf = [ax_arr[0].twinx(),ax_arr[1].twinx(),ax_arr[2].twinx()]

	ax_arr[0].plot([-99,-99],[0,0],'w.',label='[Z/H] = ')
	for nm in range(len(name_metals)):
		if name_metals[nm]=='-2.302':
			ax_arr[0].bar(-99,0,facecolor=color_metals[nm],label='-2.3 (RHB)')
		elif name_metals[nm]=='-2.301':
			ax_arr[0].bar(-99,0,facecolor=color_metals[nm],label='-2.3 (BHB)')
		elif name_metals[nm]=='-1.302':
			ax_arr[0].bar(-99,0,facecolor=color_metals[nm],label='-1.3 (RHB)')
		elif name_metals[nm]=='-1.301':
			ax_arr[0].bar(-99,0,facecolor=color_metals[nm],label='-1.3 (BHB)')
		else:
			ax_arr[0].bar(-99,0,facecolor=color_metals[nm],label=name_metals[nm])
	for i,fits in enumerate(fits_arr):

		#inc_weights = weights[fits]
		import statistics
		inc_weights = statistics.renormalise_weights(weights[fits],probs[fits])
		inc_weights = inc_weights / np.sum(inc_weights)



		prev_weight = np.zeros(len(weights[0]))
		tot_weight = np.sum(inc_weights,0)
		print tot_weight

		for mi,m in enumerate(sort_unique_metals):

			
			for w in range(len(weights[0])):


				if metal[w] == m:
					ax_arr[i].bar(age[w],tot_weight[w],width=min_width,facecolor=color_metals[mi],bottom=prev_weight[w],align='center')
					prev_weight[w] = tot_weight[w]+prev_weight[w]
		# Make a CDF:

		pdf_bar = np.zeros(len(sort_unique_ages))
		cdf_bar = np.zeros(len(sort_unique_ages))
		for a in range(len(sort_unique_ages)):
			pdf_bar[a] = np.sum(prev_weight[np.where(age==sort_unique_ages[a])])
			cdf_bar[a] = np.sum(pdf_bar[:a])
		ax_cdf[i].plot([np.min(sort_unique_ages)-99,np.min(sort_unique_ages)]+sort_unique_ages.tolist()+[np.max(sort_unique_ages),np.max(sort_unique_ages)+99],\
								[0,0]+cdf_bar.tolist()+[1.0,1.0],'--',linewidth=4.0,color='springgreen')


	# Fine-tune figure; make subplots close to each other and hide x ticks for
	# all but bottom plot.c
	ax_arr[1].set_ylabel("SSP Weight")
	ax_arr[-1].set_xlabel("Age / log(yr)")
	ax_cdf[1].set_ylabel("CDF of SSP weight")
	handles, labels = ax_arr[0].get_legend_handles_labels()

	ax_arr[0].legend(handles,labels,shadow=True,ncol=int(round(len(labels)/2)),bbox_to_anchor=(0.01, 0.95, 0.98, .105), loc=3,\
					mode="expand", borderaxespad=0.,handlelength=0.6,handleheight=1.2,handletextpad=0.4,\
					framealpha=1.0)

	ax_arr[-1].set_xlim(np.min(age[np.where(prev_weight>0.001)])-min_width,np.max(age[np.where(prev_weight>0.001)])+min_width)
	f.subplots_adjust(hspace=0)
	plt.setp([a.get_xticklabels() for a in ax_arr[:-1]], visible=False)
	plt.show()

def plot_manga():
	# Plot IFU maps
	pass
	
def files_of_light_weights(log_age, metal, weights, probs, filename=""):

	best_fit 		= np.argmax(probs)
	weights=weights
	fits_arr = [best_fit]
	for i,fits in enumerate(fits_arr):

		if i == 0:
			ii = (0,0)
			
		inc_weights = statistics.renormalise_weights(weights[fits_arr],probs[fits_arr])
		inc_weights = inc_weights / np.sum(inc_weights)
		
		prev_weight = np.zeros(len(weights[0]))
		tot_weight 	= np.sum(inc_weights,0)
		
	for i in range(len(log_age)):
				file=open('./'+str(filename)+'_age.txt', 'a')
				file.write('%12.6f\t' % log_age[i])
				file.write('\n')
				file.close()
	for i in range(len(metal)):
				file=open('./'+str(filename)+'_metal.txt', 'a')
				file.write('%12.6f\t' % metal[i])
				file.write('\n')
				file.close()
	for i in range(len(probs)):
				file=open('./'+str(filename)+'_light_probs.txt', 'a')
				file.write('%12.6f\t' % probs[i])
				file.write('\n')
				file.close()
	for i in range(len(metal)):
				file=open('./'+str(filename)+'_light_weights.txt', 'a')
				file.write('%12.6f\t' % tot_weight[i])
				file.write('\n')
				file.close()
				
def files_of_mass_weights(log_age, metal, weights, probs, filename=""):

	best_fit 		= np.argmax(probs)
	weights=weights
	fits_arr = [best_fit]
	for i,fits in enumerate(fits_arr):

		if i == 0:
			ii = (0,0)
			
		inc_weights = statistics.renormalise_weights(weights[fits_arr],probs[fits_arr])
		inc_weights = inc_weights / np.sum(inc_weights)
		
		prev_weight = np.zeros(len(weights[0]))
		tot_weight 	= np.sum(inc_weights,0)
		
	for i in range(len(log_age)):
				file=open('./'+str(filename)+'_age.txt', 'a')
				file.write('%12.6f\t' % log_age[i])
				file.write('\n')
				file.close()
	for i in range(len(metal)):
				file=open('./'+str(filename)+'_metal.txt', 'a')
				file.write('%12.6f\t' % metal[i])
				file.write('\n')
				file.close()
	for i in range(len(probs)):
				file=open('./'+str(filename)+'_mass_probs.txt', 'a')
				file.write('%12.6f\t' % probs[i])
				file.write('\n')
				file.close()
	for i in range(len(metal)):
				file=open('./'+str(filename)+'_mass_weights.txt', 'a')
				file.write('%12.6f\t' % tot_weight[i])
				file.write('\n')
				file.close()