import tempfile
import os
import distutils.core
import subprocess
import numpy as np
import shutil
import operator
import cPickle as pickle
import pylab as pl
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
                        FileTransferSpeed, FormatLabel, Percentage, \
                        ProgressBar, ReverseBar, RotatingMarker, \
                        SimpleProgress, Timer  
pl.ion()



class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

import subprocess as sub
import threading
import time as tm

class RunCmd(threading.Thread):
    def __init__(self, cmd, timeout,tempf,write_timeout):
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.timeout = timeout
        self.tempf = tempf
        self.error = 0
        self.elapsed = 0
        self.write_timeout = write_timeout

    def run(self):
        self.p = sub.Popen(self.cmd, stdout=self.tempf, stderr=self.tempf)
        self.p.wait()

    def Run(self):
        start = tm.time()
        self.start()
        self.join(self.timeout)
        self.elapsed = tm.time() - start
        
        
        if self.is_alive():
            self.p.terminate()      #use self.p.kill() if process needs a kill-9
            self.join()
            self.error = True
        if self.write_timeout:
          return self.error, self.elapsed
        else: 
          return self.error


class nautilus(object):
  """A python wrapper to nautilus"""
  def __init__(self, temproot=None, datafolder='./data', debug=False, nautilus_exec_path='/Users/gratier/code/nautilus/',timeout=60,write_timeout=False):  
    self.temproot = temproot
    self.datafolder = datafolder
    self.debug = debug
    self.timeout = timeout
    self.write_timeout = write_timeout
    self.nautilus_exec_path = nautilus_exec_path
    self.model = {}
    self.speclist = []
  
  def read_speclist(self):
    """docstring for read_speclist"""
    f = open( self.datafolder+'/gas_species.in', 'r')
    lines = f.readlines()
    for line in lines[:-3]:
      self.speclist.append(line.split()[0])
    f.close()
    f = open( self.datafolder+'/grain_species.in', 'r')
    lines = f.readlines()
    for line in lines:
      self.speclist.append(line.split()[0])
    f.close()

  def read_results(self,sigma_model):
    """docstring for read_results"""
    for fn in os.listdir('ab'):
      specname = os.path.splitext(fn)[0]
      f = open('ab/{}'.format(fn))
      lines = f.readlines()
      value = float(lines[1].split()[1])
      f.close
      self.model[specname] = (value,value/20.)
  
  def fill_results(self):
    """docstring for fill_results"""
    for specname in self.speclist:
      self.model[specname] = (np.inf,np.inf)
  
  def parameters(self):
    with cd(self.tempdir):
      parfile = 'parameters.in'
      f = open(parfile, 'r')
      lines = f.readlines()
      for line in lines:
        print line
      f.close()
  
  def write_inputs(self, time,density, gas_temperature, grain_temperature, Av, CR_rate, C_ab, N_ab, O_ab, S_ab):
    self.tempdir = tempfile.mkdtemp(prefix='tmp_', dir=self.temproot)
    distutils.dir_util.copy_tree(self.datafolder, self.tempdir)
    with cd(self.tempdir):
      parfile = 'parameters.in'
      f = open(parfile, 'r')
      lines = f.readlines()
      f.close()
      f = open(parfile, 'w')
      for line in lines:
        if len(line)>1:
          if (line.split()[0] == 'stop_time'):
            f.write('stop_time = {:.3e}\n'.format(time))
          elif (line.split()[0] == 'initial_gas_density'):
            f.write('initial_gas_density =  {:.3e}\n'.format(density))
          elif (line.split()[0] == 'initial_gas_temperature'):
            f.write('initial_gas_temperature = {:.3e}\n'.format(gas_temperature))
          elif (line.split()[0] == 'initial_dust_temperature'):
            f.write('initial_dust_temperature =  {:.3e}\n'.format(grain_temperature))       
          elif (line.split()[0] == 'initial_visual_extinction'):
            f.write('initial_visual_extinction =  {:.3e}\n'.format(Av))
          elif (line.split()[0] == 'cr_ionisation_rate'):
            f.write('cr_ionisation_rate =  {:.3e}\n'.format(CR_rate))
          elif (line.split()[0] == 'nb_outputs'):  
            f.write('nb_outputs = 1 ! Total number of outputs (used for linear or log spaced outputs)\n')
          elif (line.split()[0] == 'output_type'):  
            f.write('output_type = linear\n')
          elif (line.split()[0] == 'structure_type'):
            f.write('structure_type = 0D\n')
          elif (line.split()[0] == 'spatial_resolution'):
            f.write('spatial_resolution = 1')
          else:  
            f.write(line)
      f.close()
      
      parfile = 'abundances.in'
      f = open(parfile, 'r')
      lines = f.readlines()
      f.close()
      abundances = {}
      for line in lines[1:]:
        split =  line.split()
        specname = split[0]
        abinit = split[2]
        abinit_list = list(abinit)
        abinit_list[-4] = 'E'
        abinit = "".join(abinit_list)
        abundances[specname] = float(abinit)
        
      abundances['C+'] = C_ab
      abundances['O'] =  O_ab
      abundances['N'] =  N_ab
      abundances['S+'] = S_ab
      
      parfile = 'abundances.in'
      f = open(parfile, 'w')
      for specname in abundances.keys():
        f.write('{} = {:.3e}\n'.format(specname,abundances[specname]))
      f.close()
    
    pass
      
  def run_nautilus(self, time=1e5,density=1e5, gas_temperature=10, grain_temperature=10, Av=15, CR_rate=1.3e-17, C_ab=1.7e-4, N_ab=6.2e-5, O_ab=1.4e-4, S_ab=8e-8,write_timeout=False):
    """docstring for run_nautilus"""
    
    self.write_inputs(time,density, gas_temperature, grain_temperature, Av, CR_rate, C_ab, N_ab, O_ab, S_ab)
    with cd(self.tempdir):  
      with open(os.devnull, 'w') as tempf:
        error, elapsed = RunCmd(['{}/nautilus'.format(self.nautilus_exec_path)], self.timeout, tempf,write_timeout).Run()
        # proc = subprocess.Popen('{}/nautilus'.format(self.nautilus_exec_path), stdout=tempf, stderr=tempf)
        # proc.communicate()
        if not error:
          proc = subprocess.Popen('{}/nautilus_outputs'.format(self.nautilus_exec_path), stdout=tempf, stderr=tempf)
          proc.communicate()
          self.read_results(sigma_model=0)
        else:
          f = open('../timeout.txt','a')
          f.write("{} {} {} {} {} {} {} {} {} {} \n".format(np.log10(time),np.log10(density), gas_temperature, grain_temperature, np.log10(Av), np.log10(CR_rate), np.log10(C_ab), np.log10(N_ab), np.log10(O_ab), np.log10(S_ab)))
          f.close()
          self.fill_results()
          
    if not self.debug:
      shutil.rmtree(self.tempdir)
  
    return elapsed
    
  def sorted_model(self):
    """Sort species by increasing abundances"""
    return  sorted(self.model.items(), key=operator.itemgetter(1))
    
if __name__ == "__main__":
  N = nautilus(temproot="./",debug=True,nautilus_exec_path="/home/gratier/code/nautilus",timeout=120)
  p = [5.,5.,10., 10.,np.log10(15),np.log10(1.3e-17), np.log10(1.7e-4),np.log10(6.2e-5), np.log10(1.4e-4), np.log10(8e-8)]
  p = np.array([   5.28450033,   10.68883775,  328.40860557,   63.91882174,      1.48705375,  -16.1003161 ,   -3.95786505,   -2.0723169 ,         -5.7759699 ,   -4.88938128])
  time = pow(10,p[0])
  density = pow(10,p[1])
  gas_temperature = p[2]
  grain_temperature = p[3]
  Av = pow(10,p[4])
  CR_rate = pow(10,p[5])
  C_ab = pow(10,p[6])
  N_ab = pow(10,p[7])
  O_ab = pow(10,p[8])
  S_ab = pow(10,p[9])
  
  N.read_speclist()
  N.write_inputs(time, density, gas_temperature, grain_temperature, Av, CR_rate,C_ab, N_ab, O_ab, S_ab)
  elapsed = N.run_nautilus(time, density, gas_temperature, grain_temperature, Av, CR_rate,C_ab, N_ab, O_ab, S_ab,write_timeout=True)
  # N.parameters()
  # print N.model['CO']
  # print elapsed

  # f = open('simuldata/obs.pkl','wb')
  # pickle.dump(N.model,f)
  # f.close
  #
  # f = open('simuldata/obs.pkl','rb')
  # obs = pickle.load(f)
  # f.close
  # observed_species = ['C4H','HC3N','H2CO','NH3','l-C3H','CCS','CS','H2CCN','HC5N',
  #                     'c-C3H2','SO','CH3OH','HC7N','CH3C4H', 'CH2CHCN', 'C3S', 'C5H',
  #                      'C4H2', 'H2CS', 'C3N','H2CCO', 'CH3CN', 'c-C3H', 'HC9N', 'HCCNC',
  #                      'l-C3H2', 'HNCO', 'HCS+', 'C6H', 'C3O', 'HC3NH+', 'CH3C3N', 'CCO', 'HNCCC']
  #
  # nstep = 10
  # # steps = np.logspace(-3,4,nstep)
  # steps = np.logspace(np.log10(8e-8/1000),np.log10(1000*8e-8),nstep)
  # #
  # chi2 = []
  # #
  # pbar = ProgressBar(widgets=[Counter(), "/"+str(nstep)+" " , Percentage(), Bar(), Timer(), " ", ETA()], maxval=nstep).start()
  # cnt = 0
  # for step in steps:
  #   N.run_nautilus(S_ab=step)
  #   print N.model['CO']
  #   diff = np.array([np.log10(obs[key][0])-np.log10(N.model[key][0]) for key in observed_species])
  #   sigma = np.array([1. for key in observed_species])
  #   chi2.append(np.sum((diff/sigma)**2))
  #   cnt += 1
  #   pbar.update(cnt)
  #
  # pl.figure(1)
  # pl.clf()
  # pl.plot(steps, chi2)
  # pl.xscale('log')
  # pl.yscale('log')
  # pl.draw()
  #
  #
  #
  #
  #
