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


class nautilus(object):
  """A python wrapper to nautilus"""
  def __init__(self, temproot=None, datafolder='./data', debug=False, nautilus_exec_path='/Users/gratier/code/nautilus/'):  
    self.temproot = temproot
    self.datafolder = datafolder
    self.debug = debug
    self.nautilus_exec_path = nautilus_exec_path
    self.model = {}
  
  def read_results(self,sigma_model):
    """docstring for read_results"""
    for fn in os.listdir('ab'):
      specname = os.path.splitext(fn)[0]
      f = open('ab/{}'.format(fn))
      lines = f.readlines()
      value = float(lines[1].split()[1])
      f.close
      self.model[specname] = (value,value/2.)
  
  def run_nautilus(self, time=1e5,density=1e5, temperature=15, Av=15, CR_rate=1.3e-17, C_ab=1.7e-4, N_ab=6.2e-5, O_ab=1.4e-4, S_ab=8e-8):
    """docstring for run_nautilus"""
    self.tempdir = tempfile.mkdtemp(prefix='tmp_', dir=self.temproot)
    distutils.dir_util.copy_tree(self.datafolder, self.tempdir)
    with cd(self.tempdir):
      parfile = 'parameters.in'
      f = open(parfile, 'a')
      f.write('stop_time = {:.3e}\n'.format(time))
      f.write('initial_gas_density =  {:.3e}\n'.format(density))
      f.write('initial_gas_temperature = {:.3e}\n'.format(temperature))
      f.write('initial_visual_extinction =  {:.3e}\n'.format(Av))
      f.write('cr_ionisation_rate =  {:.3e}\n'.format(CR_rate))
      f.close()
      
      parfile = 'abundances.in'
      f = open(parfile, 'r')
      lines = f.readlines()
      f.close()
      abundances = {}
      for line in lines[1:]:
        specname, _, abinit= line.split()
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
      
      with open(os.devnull, 'w') as tempf:
        proc = subprocess.Popen('{}/nautilus'.format(self.nautilus_exec_path), stdout=tempf, stderr=tempf)
        proc.communicate()
        proc = subprocess.Popen('{}/nautilus_outputs'.format(self.nautilus_exec_path), stdout=tempf, stderr=tempf)
        proc.communicate()
      self.read_results(sigma_model=0)
    if not self.debug:
      shutil.rmtree(self.tempdir)
    
  def sorted_model(self):
    """Sort species by increasing abundances"""
    return  sorted(self.model.items(), key=operator.itemgetter(1))
    
if __name__ == "__main__":
  N = nautilus(temproot="./",debug=True,nautilus_exec_path="/home/gratier/code/nautilus")
  N.run_nautilus()
  print N.model['CO']
  f = open('simuldata/obs.pkl','wb')
  pickle.dump(N.model,f)
  f.close
  #
  # f = open('simuldata/obs.pkl','rb')
  # obs = pickle.load(f)
  # f.close
  #
  # nstep = 10
  # steps = np.logspace(1,6,nstep)
  # steps = np.linspace(0,40,nstep)
  #
  # chi2 = []
  #
  # pbar = ProgressBar(widgets=[Counter(), "/"+str(nstep)+" " , Percentage(), Bar(), Timer(), " ", ETA()], maxval=nstep).start()
  # cnt = 0
  # for step in steps:
  #   N.run_nautilus(time=step)
  #   diff = np.array([obs[key][0]-N.model[key][0] for key in obs.keys()])
  #   sigma = np.array([np.sqrt(obs[key][1]) for key in obs.keys()])
  #   chi2.append(np.sum((diff/sigma)**2))
  #   cnt += 1
  #   pbar.update(cnt)
  #
  # pl.figure(1)
  # pl.clf()
  # pl.plot(steps, chi2)
  # pl.xscale('log')
  # pl.yscale('log')
  
  
  
  
  
  