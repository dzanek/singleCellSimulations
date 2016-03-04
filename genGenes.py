import scipy.stats as ss
import scipy.special as ssp
import numpy as np
import rf, wf
import matplotlib.pylab as pl
import random as r
from drawGen import mrna_average,mrna_variance



if __name__ == '__main__':

	'''
	#make dicts from lists - {0-1, value}
	ks = dict(genSyn())
	kd = dict(genDeg())
	kon = genActiv()
	kon = dict(kon)
	'''
	params = []
	degMOD = 10
	for i in range(22000):

		nLvl = -2
	#	lookfor = r.uniform(0,1)
	#	goodKey = min(kon.keys(), key=lambda x:abs(x-lookfor))
		k_on = ss.lognorm.rvs(1.3,scale=1.3)*20**(nLvl)	# for kon and koff you can tune noise level with nLvl parameter
		k_off = ss.lognorm.rvs(0.8,loc=1,scale=2.95)*20**(nLvl)
		k_transcr = (ss.lognorm.rvs(0.95,loc=0,scale=1.6))*1 # it's set for 'per hour' but k_syn is being calculated per minute
		
		k_syn = k_transcr/(27.5*(k_on/k_off)+2.5)		# worth puting here new equation (it will give per hour rates)
		k_deg = (ss.lognorm.rvs(0.49,loc=0.01,scale=0.06))/60 # we want 'per minute' parameters so it's divided
		params.append([k_on,k_off,k_syn,k_deg,mrna_average(k_on,k_off,k_syn,k_deg),mrna_variance(k_on,k_off,k_syn,k_deg)/(mrna_average(k_on,k_off,k_syn,k_deg))**2,mrna_average(k_on,k_off,k_syn,k_deg/degMOD),mrna_variance(k_on,k_off,k_syn,k_deg/degMOD)/(mrna_average(k_on,k_off,k_syn,k_deg/degMOD))**2])
		# line above gives: parameters, expr. and var. for normal and modified k_deg
	print sum([i[4] for i in params])
	pl.cla()
	wf.writeFile('params_from_2809.pars',params) 
	#get rid of zeros for loglog plot 
	pl.loglog([i[4] for i in params if i[4] != 0 and i[5] != 0],[i[5] for i in params if i[4] != 0 and i[5] != 0],'o',label='kdeg norm')
	pl.loglog([i[6] for i in params if i[6] != 0 and i[7] != 0],[i[7] for i in params if i[6] != 0 and i[7] != 0],'o',label='kdeg / 10')
	pl.legend(loc=3)
#	pl.xlim(10**(-2),10**(5))
#	pl.ylim(0.05,10**(2))
	pl.savefig('distro.png')
	pl.cla()
	histogram = np.histogram([i[4] for i in params],100)
	pl.plot(histogram[1][1:],histogram[0],'o')
	pl.xscale('log')
        pl.yscale('log')
        pl.savefig('zipf.pdf')






