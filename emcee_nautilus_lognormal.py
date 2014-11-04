import matplotlib
matplotlib.use('Agg')


import emcee
import pylab as pl
import numpy as np
import triangle
import pynautilus 
import cPickle as pickle
import operator
import contextlib
import scipy.stats


from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
                        FileTransferSpeed, FormatLabel, Percentage, \
                        ProgressBar, ReverseBar, RotatingMarker, \
                        SimpleProgress, Timer  

@contextlib.contextmanager
def printoptions(*args, **kwargs):
    original = np.get_printoptions()
    np.set_printoptions(*args, **kwargs)
    yield 
    np.set_printoptions(**original)


def HDI(distrib,credMass=0.95449973610364158):
  n = distrib.shape[0]
  distrib_sorted = np.sort(distrib)
  exclude = n - np.floor(n * credMass)
  lowposs = distrib_sorted[0:exclude]
  uppposs = distrib_sorted[(n - exclude):n] 
  best = np.argmin(uppposs - lowposs)
  return lowposs[best], uppposs[best]
  
def confint_HDI(distrib):
  HDI_68 = HDI(distrib,credMass=scipy.special.erf(1/np.sqrt(2)))
  HDI_95 = HDI(distrib,credMass=scipy.special.erf(2/np.sqrt(2)))
  return np.array([HDI_95[0], HDI_68[0], np.median(distrib), HDI_68[1], HDI_95[1]])

def confint(x, q=[0.0228, 0.1587,.5,0.8413,0.9772]):
  return np.percentile(x, [100. * qi for qi in q])

def lnlike(p,obs):
  nobs = len(obs)
  time = pow(10,p[0])
  density = pow(10,p[1])
  temperature = p[2]
  Av = pow(10,p[3])
  CR_rate = pow(10,p[4])
  C_ab = pow(10,p[5])
  N_ab = pow(10,p[6])
  O_ab = pow(10,p[7])
  S_ab = pow(10,p[8])
  # print time, density, temperature, Av
  if (limits[0:ndim,0]<=p).all() and (limits[0:ndim,1]>=p).all():
    N.run_nautilus(time, density, temperature, Av, CR_rate,C_ab, N_ab, O_ab, S_ab)
    # N.run_nautilus(time, density, temperature, Av)
    
    diff = np.array([obs[key][0]-np.log10(N.model[key][0]) for key in observed_species])
    sigma = np.array([obs[key][1] for key in observed_species])
    lnlike = -0.5*np.sum((diff/sigma)**2)-0.5*nobs*np.log(2*np.pi)-np.sum(np.log(sigma))
  else:
    lnlike = -np.inf
  return lnlike, N.model

def lnlike_noblob(p, obs):
  return lnlike(p,obs)[0]
  
def lnprior(p):
  if (limits[0:ndim,0]<=p).all() and (limits[0:ndim,1]>=p).all():
    lnprior = scipy.stats.norm.logpdf(p[2],loc=10.,scale=5.)
  else:
    lnprior = -np.inf
  return lnprior
  
def lnprob(p,obs):
  loglike, blobs = lnlike(p,obs)
  return loglike+lnprior(p), blobs



f = open('simuldata/obs.pkl','rb')
obs = pickle.load(f)
f.close
for name in obs.keys():
  obs[name] = (np.log10(obs[name][0]), .3)

observed_species = ['C4H','HC3N','H2CO','NH3','l-C3H','CCS','CS','H2CCN','HC5N', 
                    'c-C3H2','SO','CH3OH','HC7N','CH3C4H', 'CH2CHCN', 'C3S', 'C5H',
                     'C4H2', 'H2CS', 'C3N','H2CCO', 'CH3CN', 'c-C3H', 'HC9N', 'HCCNC',
                     'l-C3H2', 'HNCO', 'HCS+', 'C6H', 'C3O', 'HC3NH+', 'CH3C3N', 'CCO', 'HNCCC']


# f = open('tmc1_obs_bayes.txt','r')
# lines = f.readlines()
# f.close()
# obs = {}
# observed_species = []
# for line in lines:
#   specname, column_density, sigma = line.split()
#   obs[specname] = (float(column_density)-22, float(sigma))
#   observed_species.append(specname)
  

f = open('tmc1_obs_ohishi.txt','r')
lines = f.readlines()
f.close()
obs_ohishi = {}
for line in lines:
  specname, column_density= line.split()
  obs_ohishi[specname] = (float(column_density)-22.3, 0.3)

obs = obs_ohishi


N = pynautilus.nautilus(temproot="./",debug=False,nautilus_exec_path="/home/gratier/code/nautilus")
nobs = len(obs)
labels = np.array(["$\log t (yr)$", "$\log n_H (cm^{-3}$)", "$T (K)$", "$\log Av (mag)$", "$\log \zeta (s^{-1})$","$\log C/H $","$\log N/H $","$\log O/H $","$\log S/H $"])
pinit = [5.,5.,15.,np.log10(15)]
limits = np.array([[1.,8.],
                   [1.,10.],
                   [5.,100.],
                   [-2.,2.],
                   [-18.,-15.],
                   [-10,-2],
                   [-10,-2],
                   [-10,-2],
                   [-10,-2]])
                       
ntemps = 3
ndim = 9
nwalkers = 20
nburn = 2000
niter = 50
nthreads = 24
nthin = min(nburn,10)
prefix = "lognorm"

restart_from_previous_run = True 
if restart_from_previous_run:
  f = open("chains/"+prefix+"_chain.dat", "rb")
  p0 = pickle.load(f)
  f.close()
else:
  p0 = np.random.uniform(low=limits.T[0,0:ndim],high=limits.T[1,0:ndim], size=(ntemps, nwalkers, ndim))

PTsampler = emcee.PTSampler(ntemps, nwalkers, ndim, logl=lnlike_noblob, logp=lnprior, loglargs=[obs],threads = nthreads)
pos = np.zeros((nwalkers,ndim))
lnprobval = np.zeros((nwalkers))

if nburn!=0:
  cnt = 0
  pbar = ProgressBar(widgets=[Counter(), "/"+str(nburn)+" " , Percentage(), Bar(), Timer(), " ", ETA()], maxval=nburn).start()
  for result in PTsampler.sample(p0, iterations=nburn, storechain=True, thin=nburn/nthin):
    cnt += 1
    pos = result[0]
    lnprobval = result[1]
    f = open("chains/"+prefix+"_chain.dat", "wb")
    pickle.dump(pos,f)
    f.close()
    pbar.update(cnt)

  pl.figure(2)
  pl.clf()
  orderedlnprobabilityidx = np.argsort(lnprobval[0])
  ordered = -lnprobval[0][orderedlnprobabilityidx[::-1]]
  lht = (ordered[1:]-ordered[:-1])
  rht = 5*(ordered - ordered[0])/(np.arange(nwalkers)+1)
  cut = np.argmax(lht > rht[1:])
  pl.subplot(221)
  pl.plot(ordered,np.arange(nwalkers))
  if cut:
    pl.hlines(cut,0,ordered[np.isfinite(ordered)][-1])
  pl.subplot(222)
  pl.plot(ordered,np.arange(nwalkers))
  if cut:
    pl.hlines(cut,0,ordered[np.isfinite(ordered)][-1])
    pl.ylim((0,cut))
    pl.xlim((ordered[0],ordered[cut]))
  pl.subplot(223)
  pl.plot(lht,np.arange(nwalkers)[1:])
  pl.plot(rht,np.arange(nwalkers))
  if cut:
    pl.hlines(cut,0,np.min((lht[np.isfinite(lht)].min(),lht[np.isfinite(lht)].min())))
  pl.subplot(224)
  pl.plot(lht,np.arange(nwalkers)[1:])
  pl.plot(rht,np.arange(nwalkers))
  if cut:
    pl.ylim((0,cut+10))
    pl.xlim((0,np.min((lht[np.isfinite(lht)][0:cut+10].min(),lht[np.isfinite(lht)][0:cut+10].min()))))
    pl.hlines(cut,0,np.min((lht[np.isfinite(lht)].min(),lht[np.isfinite(lht)].min())))
  pl.draw()
  pl.savefig("ha/"+prefix+"_lnprob.jpg")
    
  Z = np.zeros((nthin,3,np.min([30,ndim])))
  cnt = 0
  for k in range(np.min([30,ndim])):
    for i in range(nthin):
      Z[i,:,k] = np.array(confint(PTsampler.chain[0,:,i,k])[1:4])
  
  pl.figure(3)
  pl.clf()    
  for i in range(np.min([30,ndim])):
    pl.subplot(3,3,i+1)
    pl.plot(Z[:,0,i],'r-')
    pl.plot(Z[:,1,i],'b-')
    pl.plot(Z[:,2,i],'r-')
    pl.ylim(limits[i,:])
  pl.draw()
  pl.savefig("ha/"+prefix+"_burnin_convergence.jpg")
  pl.close(2)
  
  for itemp in range(ntemps):
    fig2 = triangle.corner(PTsampler.chain[itemp,:,-1,:].reshape(-1,ndim)[:,0:12], extents=limits[0:12], labels=labels[0:12],quantiles=[0.02275,0.1586,0.5,0.8413,0.977249],verbose=False,plot_datapoints=False)
    pl.savefig("ha/"+prefix+"_corner_"+str(itemp)+"_final.jpg")
  
  print("Mean acceptance fraction:", np.mean(PTsampler.acceptance_fraction))
  
  print("Estimate of the evidence:", PTsampler.thermodynamic_integration_log_evidence(fburnin = 0.1))
  
if niter!=0:
  pos1 = pos[0]
  sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[obs],threads = nthreads)
  cnt = 0
  pbar = ProgressBar(widgets=[Counter(), "/"+str(niter)+" " , Percentage(), Bar(), Timer(), " ", ETA()], maxval=niter).start()
  for result in sampler.sample(pos1, iterations=niter, storechain=True):
    cnt += 1
    pbar.update(cnt)

  print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))
  percentiles = np.array([2.27501319,  15.86552539,  50.        ,  84.13447461,  97.72498681])
  
  flat = sampler.chain[:,:,:].reshape(-1,ndim)
  triangle.corner(flat[:,0:12], extents=None, labels=labels[0:12],bins=15,truths=None,show_titles=True,quantiles=percentiles[1:-1]/100.)
  pl.savefig("ha/"+prefix+"_corner.jpg")
  pl.close()
  
  with printoptions(precision=3, suppress=True):
    output = np.zeros((ndim,5))
    print "confint"
    for i in range(ndim):
      output[i,:] = np.array(confint(flat[:,i]))
      print "{:3} {:10} {} {}".format(i, labels[i], output[i,:], flat[:,i].mean())
      #print i, labels[i], confint_HDI(flat[:,i]), flat[:,i].mean()
      
      
  flatblobs = np.array(sampler.blobs).reshape(-1)
  N.run_nautilus()

  dict_blob = {}
  for specname in N.model.keys():
    tmp = []
    for blob in flatblobs:
      tmp.append(np.log10(blob[specname][0]))
    dict_blob[specname] = np.array(np.percentile(tmp,q=percentiles,axis=0))+0.3

  
  sorted_ab = sorted(dict_blob.iteritems(), key=lambda (k,v): operator.itemgetter(2)(v)) 

  pl.figure(4,figsize=(1,10))
  pl.clf()
  pos = np.arange(len(observed_species))+.5
  spec_labels = []
  for i, specname in enumerate(observed_species[::-1]):
    pl.barh(bottom = pos[i],align='center', width = dict_blob[specname][-1]-dict_blob[specname][0],left= dict_blob[specname][0],edgecolor='none',color='k', alpha=0.2)
    pl.barh(bottom = pos[i],align='center', width = dict_blob[specname][-2]-dict_blob[specname][1],left= dict_blob[specname][1],edgecolor='none',color='k', alpha=0.4)
    pl.errorbar(obs[specname][0]+0.3, pos[i], xerr=obs[specname][1],color='r')
    if obs_ohishi[specname][0]<0:
      pl.scatter(obs_ohishi[specname][0]+0.3, pos[i],marker='|',linewidth=3,s=100)
    spec_labels.append(specname)

  pl.yticks(pos, spec_labels)  
  pl.grid(b='on',which='both',axis='both')
  pl.xlabel('Abundance relative to $H_2$')
  pl.savefig('ha/model.png')


  fig = pl.figure(5)
  pl.clf()
  pos = 2*np.arange(len(sorted_ab))+0.5
  spec_labels = []
  i = 0
  for value in sorted_ab:
    specname, confint = value
    if confint[0] > -20:
      pl.barh(bottom = pos[i],align='center', width = confint[-1]-confint[0],left= confint[0],edgecolor='none',color='k', alpha=0.2)
      pl.barh(bottom = pos[i],align='center', width = confint[-2]-confint[1],left= confint[1],edgecolor='none',color='k', alpha=0.4)
      spec_labels.append(specname)
      i += 1


  from matplotlib.ticker import MultipleLocator, FormatStrFormatter

  majorLocator   = MultipleLocator(10)
  majorFormatter = FormatStrFormatter('%d')
  minorLocator   = MultipleLocator(2)

  pl.yticks(pos[0:i], spec_labels)  
  pl.xlabel('Abundance relative to $H_2$')
  pl.gca().xaxis.set_major_locator(majorLocator)
  pl.gca().xaxis.set_major_formatter(majorFormatter)
  pl.gca().xaxis.set_minor_locator(minorLocator)
  pl.gca().xaxis.set_tick_params(labeltop='on')
  pl.gca().yaxis.set_tick_params(labelright='on')
  
  pl.grid(b='on',which='both',axis='both')
  fig.set_size_inches(50,100)
  pl.ylim((pos[0],pos[-1]))
  pl.savefig('ha/full_model.png')
  