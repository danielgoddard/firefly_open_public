#parameters_obtain.py

def parameters_obtain(parameters_file='./parameters.dat'):

	"""
	Returns the a dictionary of options that are easily python-accessible
	from the options file parameters.py. Unspecified (commented-out) options are set to defaults.
	"""

	#######################
	### DEFAULT OPTIONS: DO NOT EDIT ###

	metal_used 			= ['z001','z002','z004','z0001.bhb','z0001.rhb','z10m4.bhb','z10m4.rhb']
	number_metal_used 	= [0.5,1.0,2.0,0.05,0.05,0.005,0.005]

	age_used 			= ['3M','3_5eM','4M','4_5eM','5M','5_5eM','6M','6_5eM',\
							'7M','7_5eM','8M','8_5eM','9M','9_5eM',\
							'10M','15M','20M','25M','30M','35M','40M','45M','50M','55M',\
							'60M','65M','70M','75M','80M','85M','90M','95M',\
							'100M','200M','300M','400M','500M','600M','700M','800M','900M',\
							'1G','1G_500M','2G','3G','4G','5G','6G','7G','8G','9G',\
							'10G','11G','12G','13G','14G','15G']
	number_age_used 	= [3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,\
							10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,\
							100,200,300,400,500,600,700,800,900,\
							1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,\
							10000,11000,12000,13000,14000,15000]

	model_libs 			= ['MILES','STELIB','ELODIE','MILES_revisedIRslope','MILES_UVextended']
	imfs 				= ['kr','ss','cha']

	default_options 	= {	'data_type' : 'unspecified',\
							'file_in' 	: 'unspecified',\
							'models'	: 'm11',\
							'age_used'  : age_used,\
							'metal_used': metal_used,\
							'number_age_used':number_age_used,\
							'number_metal_used':number_metal_used,\
							'model_libs' : model_libs,\
							'imfs'		: imfs,\
							'min_wavelength': 0,\
							'max_wavelength':0,\
							'hpf_mode':'on'}


	### ^ DO NOT EDIT ^ ###
	########################

	import parameters

	# Adapted from spouk's answer on Stack Overflow:
	# http://stackoverflow.com/questions/9759820/how-to-get-a-list-of-variables-in-specific-python-module
	options = {key: value for key, value in parameters.__dict__.iteritems() if not (key.startswith('__') or key.startswith('_'))}
	
	# print "------------------------------------------------------"
	# print "List of parameters using to fit:"
	# print "------------------------------------------------------"
	for key,value in default_options.iteritems():
		# If parameter not specified in file, set to default:
		if not options.has_key(key):
			options[key] = value
		#print key+' = '+str(value)
	print "------------------------------------------------------"

	if options['file_in'] == 'unspecified' or options['data_type'] == 'unspecified':
		print "You have not specified a file!!"
		print "help..."

	print options
	return options
	 