import pylab as pl
import numpy as np
import scipy.stats as stats

pl.ion()


normpdf = lambda x,p: 1/np.sqrt(2*np.pi*p[1]**2)*np.exp(-((x-p[0])**2/(2*p[1]**2)))
x = np.logspace(-1,3,10000)
y = normpdf(x,[0,1])
y1 = stats.lognorm.pdf(x,0.5,10)

pl.figure(1)
pl.clf()
pl.plot(x,y)
pl.plot(x,y1)

pl.xscale('log')
pl.yscale('linear')