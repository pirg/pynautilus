import matplotlib
matplotlib.use('Agg')

import emcee
import pylab as pl
import numpy as np
import triangle
import pynautilus 
import cPickle as pickle
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
                        FileTransferSpeed, FormatLabel, Percentage, \
                        ProgressBar, ReverseBar, RotatingMarker, \
                        SimpleProgress, Timer  


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
  # print time, density, temperature, Av
  if (limits[:,0]<=p).all() and (limits[:,1]>=p).all():
    model = N.run_nautilus(time, density, temperature, Av)
    diff = np.array([obs[key][0]-N.model[key][0] for key in obs.keys()])
    sigma = np.array([np.sqrt(obs[key][1]) for key in obs.keys()])
    lnlike = -0.5*np.sum((diff/sigma)**2)#-0.5*nobs*np.log(2*np.pi)-np.sum(np.log(sigma))
  else:
    lnlike = -np.inf
  return lnlike

  
def lnprior(p):
  if (limits[:,0]<=p).all() and (limits[:,1]>=p).all():
    lnprior = 0
  else:
    lnprior = -np.inf
  return lnprior
  
def lnprob(p,obs):
  loglike = lnlike(p,obs)
  return loglike+lnprior(p)



f = open('simuldata/obs.pkl','rb')
obs = pickle.load(f)
f.close

N = pynautilus.nautilus(temproot="./",debug=False,nautilus_exec_path="/home/gratier/code/nautilus")
nobs = len(obs)
labels = np.array(["time", "density", "temperature", "Av"])
pinit = [5.,5.,15.,np.log10(15)]
limits = np.array([[1.,8.],
                   [1.,10.],
                   [5.,100.],
                   [-2,2]])
                       
ntemps = 5
ndim = 4
nwalkers = 200
nburn = 1
niter = 1
nthreads = 16
nthin = min(nburn,10)
prefix = ""

restart_from_previous_run = True
if restart_from_previous_run:
  f = open("chains/"+prefix+"_chain.dat", "rb")
  p0 = pickle.load(f)
  f.close()
else:
  p0 = np.random.uniform(low=limits.T[0],high=limits.T[1], size=(ntemps, nwalkers, ndim))

PTsampler = emcee.PTSampler(ntemps, nwalkers, ndim, logl=lnlike, logp=lnprior, loglargs=[obs],threads = nthreads)
pos = np.zeros((nwalkers,ndim))
lnprobval = np.zeros((nwalkers))

# nstep = 20
# steps = np.linspace(1,10,nstep)
# # steps = np.linspace(0,40,nstep)
#
# chi2 = []
#
# pbar = ProgressBar(widgets=[Counter(), "/"+str(nstep)+" " , Percentage(), Bar(), Timer(), " ", ETA()], maxval=nstep).start()
# cnt = 0
# for step in steps:
#   # N.run_nautilus(time=pow(10,step))
#   # diff = np.array([obs[key][0]-N.model[key][0] for key in obs.keys()])
#   # sigma = np.array([np.sqrt(obs[key][1]) for key in obs.keys()])
#   # chi2.append(np.sum((diff/sigma)**2))
#
#   chi2.append(-lnlike([step,5.,15.,np.log10(15)],obs))
#   cnt += 1
#   pbar.update(cnt)
#
#
# pl.figure(1)
# pl.clf()
# pl.plot(steps, chi2)
# pl.xscale('linear')
# pl.yscale('linear')
# pl.ylim([0,3])




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
    pl.subplot(2,2,i+1)
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
  
  flat = sampler.chain[:,:,:].reshape(-1,ndim)
  triangle.corner(flat[:,0:12], extents=limits[0:12], labels=labels[0:12],bins=20,truths=pinit)
  pl.savefig("ha/"+prefix+"_corner.jpg")
  pl.close()
  
  # with printoptions(precision=3, suppress=True):
  #   output = np.zeros((ndim,5))
  #   print "confint"
  #   for i in range(ndim):
  #     output[i,:] = np.array(confint(flat[:,i]))
  #     print "{:3} {:10} {} {}".format(i, labels[i], output[i,:], flat[:,i].mean())
  #     #print i, labels[i], confint_HDI(flat[:,i]), flat[:,i].mean()

  