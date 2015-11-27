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
from matplotlib.colors import LogNorm

def main():

	#data=open('/Users/Daniel/test_output.txt', 'r')
	#light_age=[]
	#light_metal=[]
	#light_probs=[]
	#light_weights=[]

#	for line in data.readlines():
#		ln=line.split()
#		age=float(ln[0])
#		metal=float(ln[1])
#		probs=float(ln[2])
#		weights=float(ln[3])
#		light_age.append(age)
#		light_metal.append(metal)
	#	light_probs.append(probs)
	#	light_weights.append(weights)
	#light_age=np.array(light_age)
	#light_metal=np.array(light_metal)
	#light_probs=np.array(light_probs)
	#light_weights=np.array(light_weights)
	
	data=open('/Users/Daniel/firefly_original/mass_output_age.txt', 'r')
	light_age=[]
	for line in data.readlines():
		ln=line.split()
		age=float(ln[0])
		light_age.append(age)
	light_age=np.array(light_age)
	
	data=open('/Users/Daniel/firefly_original/mass_output_metal.txt', 'r')
	light_metal=[]
	for line in data.readlines():
		ln=line.split()
		metal=float(ln[0])
		light_metal.append(metal)
	light_metal=np.array(light_metal)
	
	data=open('/Users/Daniel/firefly_original/mass_probs.txt', 'r')
	light_probs=[]
	for line in data.readlines():
		ln=line.split()
		probs=float(ln[0])
		light_probs.append(probs)
	light_probs=np.array(light_probs)
	
	data=open('/Users/Daniel/firefly_original/mass_tot_weight.txt', 'r')
	light_weights=[]
	for line in data.readlines():
		ln=line.split()
		weights=float(ln[0])
		light_weights.append(weights)
	light_weights=np.array(light_weights)
	#print np.shape(light_weights)
	#print np.shape(light_probs)
	#print np.shape(light_metal)
	#print np.shape(light_age)
	
	
	
	light_weights=light_weights[0:401498]
	light_probs=light_probs
	light_metal=light_metal[0:401498]
	light_age=light_age[0:401498]


	''' define variables '''
	age = light_age+9.0
	metal = np.asarray(light_metal)
	probs=light_probs
	#weights=light_weights
	
	def diff_array(array):
		bisected_array = np.zeros(len(array)-1)
		for ai in range(len(bisected_array)):
			bisected_array[ai] = (array[ai+1]-array[ai])
		return bisected_array

	min_width = np.min(np.absolute((diff_array(age))))
	

	''' Plots the single contour Dan version.'''
	best_fit 		= np.argmax(probs)
	fits_arr=[best_fit]
	label_arr = ["DR7 Sub-Sample"]
	col_cum = ['black']
	sort_ages 	= np.sort(age)
	sort_age_index = np.argsort(age)
	sort_unique_ages = np.sort(np.unique(age))
	sort_metals = np.sort(metal)
	sort_unique_metals = np.sort(np.unique(metal))
	name_metals = [str(x) for x in sort_unique_metals]

	
	def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
		new_cmap = colors.LinearSegmentedColormap.from_list(
		'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
		cmap(np.linspace(minval, maxval, n)))
		return new_cmap

	cmap=truncate_colormap(get_cmap('hot_r'),minval=0.15,maxval=1.0,n=1000)

	cmap.set_under('w')
	bounds_interp_age = (9.0,10.19)
	bounds_interp_metal = (-2.5,0.90)
	
	
	f=plt.figure(1)
	plt.tick_params(labelsize=14)
	plt.plot([0,0],[0,0],'-',color=col_cum[0],markersize=1,label=label_arr[0],linewidth=3.0)
	tot_weight=light_weights

	index_allow = np.where( (age>0) & (metal>-90) )
	npts 	= 300
	extent_array = [bounds_interp_age[0],bounds_interp_age[1],bounds_interp_metal[0],bounds_interp_metal[1]]

	xage 	= np.linspace(bounds_interp_age[0],bounds_interp_age[1],npts)
	ymetal 	= np.linspace(bounds_interp_metal[0],bounds_interp_metal[1],npts)
	plt.xlim(bounds_interp_age[0],bounds_interp_age[1])
	plt.ylim(bounds_interp_metal[0],bounds_interp_metal[1])
	# grid the data.
	zi = griddata((age, metal), tot_weight, (xage[None,:], ymetal[:,None]), method='linear')
	zi[np.isnan(zi)|np.isinf(zi)]=0.00001
	zs =scipy.ndimage.gaussian_filter(zi,sigma=5,order=0)
	
	CS = plt.contour(xage,10**ymetal,zs/np.max(zs),[0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],\
								linewidths=0.5,colors='k',vmin=0.01)
	CS = plt.contourf(xage,10**ymetal,zs/np.max(zs),100,cmap=cmap,vmin=0.01,vmax=1.0)
	plt.legend(loc="upper left",prop={'size':18},shadow=False)

	plt.ylabel(r"$[Z/H]$",fontsize=24)
	plt.xlabel("Age($\log(yr)$)",fontsize=20)
	plt.ylim(0.001,2.0)
	plt.xlim(9.4,10.2)
	plt.tight_layout()
	cb = f.colorbar(CS,cmap=cmap)
	cb.set_ticklabels(['0.01','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'])
	cb.set_label("Relative SSP weight",fontsize=20,labelpad=5)
	cb.ax.tick_params(labelsize=18)
	
	plt.show()
	

main()
	