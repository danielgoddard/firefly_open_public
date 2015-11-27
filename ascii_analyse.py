#ascii_analyse.py

import os
import numpy as np
import matplotlib.pyplot as plt

def ascii_analyse():

	file_dir = 'fits/moresco/'
	model_select = 'STELIB_BC03/cha/'
	restrict_age = False

	if restrict_age:
		len_str = 10
		sub_len = 0
		fit_dir = file_dir+"moresco_allages/"+model_select
		age_str = "_restrict_age"
	else:
		len_str = 9
		sub_len = 1
		fit_dir = file_dir+"moresco_restrictage/"+model_select
		age_str = "_all_ages"



	redshift,vdisp,light_age,light_metal,mass_age,mass_metal,mass = \
	[],[],[],[],[],[],[]

	for f in os.listdir(fit_dir):
		print f
		if f[0]!='a':
			continue
		print f

		file_split = f.split('_')


		print file_split
		vdisp_str = file_split[1]
		if vdisp_str == 'mm1':
			vdisp.append(100.0)
		elif vdisp_str == 'mm2':
			vdisp.append(175.0)
		elif vdisp_str == 'mm3':
			vdisp.append(225.0)
		elif vdisp_str == 'mm4':
			vdisp.append(275.0)
		elif vdisp_str == 'mm5':
			vdisp.append(350.0)
		redshift.append(float(file_split[3][1:]))
		fits = np.loadtxt(fit_dir+f, skiprows=1, unpack=True)
		light_age.append(fits[0])
		light_metal.append(fits[1])
		mass_age.append(fits[2])
		mass_metal.append(fits[3])
		mass.append(fits[4])


	plt.plot(redshift,light_age,'bo')
	plt.xlabel("Redshift",fontsize=22)
	plt.ylabel("Light-weighted age / Gyr",fontsize=22)
	out_plot_string = 'plots/moresco_redshift_lightage'+age_str+'.eps'
	plt.savefig(out_plot_string,format='eps',transparent=True)
	plt.close()

	plt.plot(redshift,light_metal,'go')
	plt.xlabel("Redshift",fontsize=22)
	plt.ylabel("Light-weighted [Z/H]",fontsize=22)
	out_plot_string = 'plots/moresco_redshift_lightmetal'+age_str+'.eps'
	plt.savefig(out_plot_string,format='eps',transparent=True)
	plt.close()


	plt.plot(vdisp,light_age,'bo')
	plt.xlabel("Velocity dispersion / km/s",fontsize=22)
	plt.ylabel("Light-weighted age / Gyr",fontsize=22)
	out_plot_string = 'plots/moresco_vdisp_lightage'+age_str+'.eps'
	plt.savefig(out_plot_string,format='eps',transparent=True)
	plt.close()
	
	plt.plot(vdisp,light_metal,'go')
	plt.xlabel("Velocity dispersion / km/s",fontsize=22)
	plt.ylabel("Light-weighted [Z/H]",fontsize=22)
	out_plot_string = 'plots/moresco_vdisp_lightmetal'+age_str+'.eps'
	plt.savefig(out_plot_string,format='eps',transparent=True)
	plt.close()


	plt.plot(redshift,mass_age,'bo')
	plt.xlabel("Redshift",fontsize=22)
	plt.ylabel("Mass-weighted age / Gyr",fontsize=22)
	out_plot_string = 'plots/moresco_redshift_massage'+age_str+'.eps'
	plt.savefig(out_plot_string,format='eps',transparent=True)
	plt.close()
	plt.plot(redshift,mass_metal,'go')
	plt.xlabel("Redshift",fontsize=22)
	plt.ylabel("Mass-weighted [Z/H]",fontsize=22)
	out_plot_string = 'plots/moresco_redshift_massmetal'+age_str+'.eps'
	plt.savefig(out_plot_string,format='eps',transparent=True)
	plt.close()

	plt.plot(vdisp,mass_age,'bo')
	plt.xlabel("Velocity dispersion / km/s",fontsize=22)
	plt.ylabel("Mass-weighted age / Gyr",fontsize=22)
	out_plot_string = 'plots/moresco_vdisp_massage'+age_str+'.eps'
	plt.savefig(out_plot_string,format='eps',transparent=True)
	plt.close()
	plt.plot(vdisp,mass_metal,'go')
	plt.xlabel("Velocity dispersion / km/s",fontsize=22)
	plt.ylabel("Mass-weighted [Z/H]",fontsize=22)
	out_plot_string = 'plots/moresco_vdisp_massmetal'+age_str+'.eps'
	plt.savefig(out_plot_string,format='eps',transparent=True)
	plt.close()










