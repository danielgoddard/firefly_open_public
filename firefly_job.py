import numpy as np
from firefly_single import *
import glob
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

from IPython.core.debugger import Tracer

def firefly_job(job_begin,job_end):

  """
  This is a wrapper function for running firefly_single
  on multiple data SEDs. For more information on this process,
  see firefly_single.py or the README.

  Given a directory, it analyses all data at the locations
  between [job_begin] (inclusive) percent and [job_end] (exclusive) 
  percent of all the files in the directory.
  e.g. if 1000 files are given and the command is:
  firefly_job([file_dir],10,11)
  FIREFLY will analyses file numbers 100 to 109 in that directory.


  Input:
  - data_dir:   directory of data SEDs to be analysed
  - job_begin:  starting percentage of files to analyse (inclusive)
  - job_end:    finishing percentage of files to analyse (exclusive)


  """

  parameters  = parameters_obtain('parameters.py')

  input_dir   = parameters['file_dir']
  file_list   = os.listdir(input_dir)#glob.glob(input_dir+'*.dat')

  num_files   = len(file_list)

  start_job   = num_files / 100.0 * job_begin
  end_job   = num_files / 100.0 * job_end


  if parameters['observation_type'] == 'ifu':
    count_ifu = 0
    file_ifu,bin_ifu = [],[]
    for f in range(num_files):
      parameters['file_in']     = file_list[f]

      if file_list[f][0]=='.':
        continue
      #Need to loop over bins within the ifu datacube as well as file.
      header      = pyfits.open(parameters['file_dir']+parameters['file_in'],ignore_missing_end=True)
      flux_all    = header['FLUX'].data
      maxshape    = np.shape(flux_all)[1]
      
      #len_ifu += maxshape
      for l in range(maxshape):
        file_ifu.append(parameters['file_dir']+file_list[f])
        bin_ifu.append(l)
        count_ifu += 1

    start_job   = count_ifu / 100.0 * job_begin
    end_job   = count_ifu / 100.0 * job_end
    #Tracer()()
    for ff in range(count_ifu):

      if ff < end_job and ff >= start_job:

        parameters['file_in']     = file_ifu[ff]


        parameters['bin_number']  = bin_ifu[ff]

        print "------------------------------------------------------"
        print "List of parameters using to fit:"
        print "------------------------------------------------------"
        for key,value in parameters.iteritems():
          print key+' = '+str(value)
        print "------------------------------------------------------"
        firefly_single(parameters)

        
  else:
    for f in range(num_files):

      if f < end_job and f >= start_job:

        print input_dir
        print file_list
        if file_list[f][0] == '.':
          continue

        parameters['file_in']     = input_dir+file_list[f]
        parameters['output_file']   = parameters['output_dir_prefix']+parameters['file_in'].split('/')[-1]



        print "------------------------------------------------------"
        print "List of parameters using to fit:"
        print "------------------------------------------------------"
        for key,value in parameters.iteritems():
          print key+' = '+str(value)
        print "------------------------------------------------------"
        firefly_single(parameters)