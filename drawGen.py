import numpy as np
import scipy.special
import math
import scipy.integrate
import pandas as pd
from scipy import interpolate
from scipy.optimize import  minimize_scalar
from sys import argv as av


def readFile(fname):
	plik = open(fname,'r').readlines()
	plik = [i.split() for i in plik]
	out = []
	for i in plik:
		item = []
		for k in i:
			try: item.append(float(k))
			except: item.append(str(k))
		if len(item) > 0:
			out.append(item)


#	plik = [[float(k) for k in i] for i in plik]
	return out



def mrna_average(k0=1.0,k1=1.0,v0=8.0,d0=0.1):#{{{
    return k0*v0/((k0+k1)*d0)#}}}

def mrna_variance(k0=1.0,k1=1.0,v0=8.0,d0=0.1):#{{{
    va=(k0*v0/(d0*(k0+k1)))*(1+(k1*v0/((k0+k1)**2+d0*(k0+k1))))
    return va#}}}

def mrna_probability(m=np.array(1.0),l=1.0,g=1.0,mu=8.0,d=0.1):#{{{
#taken from supplementary of Raj Plos biol 2006
    result=[[],[]]
    for i in m:
        try:
            #pdb.set_trace()
            f1=scipy.special.hyp1f1(l/d + i, l/d + g/d +i , -1*mu/d)
            t1=math.gamma(l/d+i)/math.gamma(i+1)
            t2=(mu/d)**i/math.gamma(l/d+g/d+i)
            t3=math.gamma(l/d + g/d)/math.gamma(l/d)
            p=t1*t2*t3*f1
            if p != float('NaN')  and p != float('Inf') and p > 0:
                result[0].append(i)
                result[1].append(p)
        except OverflowError:
            pass
    return np.array(result)#}}}

def spline_fit_function(x=0,y=0,inv_cumulative=0,der=0):#{{{
    return abs(interpolate.splev(x,inv_cumulative,der) -y)#}}}

def bin_data(coordinate,data, nr_bins=100):#{{{
  "bins some data in a given direction (coordinate), e.g. x"
  bins=np.linspace(coordinate.min(), coordinate.max(), nr_bins+1)
  digitized = np.digitize(coordinate, bins)
  #pdb.set_trace()
  bin_means = [data[digitized == i].mean() for i in range(1, (len(bins)-1))]
  #bin_var = [data[digitized == i].var() for i in range(1, (len(bins)))]
  bin_mean_pos=[]
  for i in range(len(bins)-1):
      bin_mean_pos.append(np.average((bins[i+1], bins[i])))

  return np.array(bin_mean_pos),np.array(bin_means)
#}}}

def mrna_sample_probability(l=1.0,g=1.0,mu=8.0,d=0.1, the_sample=1000):#{{{
#Sample the probability using numerical inverse cdf aproach

    #set the scale of the distribution
    average=mrna_average(k0=l,k1=g,v0=mu,d0=d)
    variance=mrna_variance(k0=l,k1=g,v0=mu,d0=d)
    print "Average: ", average, "Variance: ", variance
    #0.001 is arbitrary
    x=np.arange(max(0,average-variance), average+variance, 0.01)

    #calculate probability didtribution
    distribution=mrna_probability(x,l,g,mu,d)

    #calculate cumulative function, returns one value less, i.e. without x[0]
    cumulative=scipy.integrate.cumtrapz(distribution[1], distribution[0],initial=0)

    #cumulative has problems with high values on zero
    cumulative+=1-cumulative[-1]
    cumulative=pd.DataFrame(np.array([distribution[0],cumulative]).T)
    cumulative=cumulative[cumulative[1]>0]
    cumulative=cumulative[:np.where(cumulative[1]==1.0)[0][0]]
    if cumulative.iloc[0][0] == 0.0:
        the_zero=cumulative[0][1]
        cumulative=cumulative[1:]
    #log space ?
    #cumulative[0]+=1
    #cumulative=np.log(cumulative)
    binned_cum=bin_data(cumulative[0] , cumulative[1] , 200)
    #print "Range: ", np.exp(cumulative[0].min())-1,",", np.exp(cumulative[0].max())-1
    #print "Probability integral: ", np.exp(cumulative[1][-1:].values)
 #   print "Range: ", cumulative[0].min(),",", cumulative[0].max()
#    print "Probability integral: ", cumulative[1][-1:].values

    #fit inverse cumulative function
    #check binned , why one is one less

    inv_cumulative =interpolate.splrep(binned_cum[0][:-1],binned_cum[1],s=0) #is s = 0 no smoothing?

    #sample the values of the distribution using the inverse cumulative and relation to uniform distribution
    #random is [0,1), should it be [0,1]?
    #random_values=[]
    #for i in range(the_sample):
    #    random_values.append(random.uniform(cumulative[1].min(), cumulative[1].max()))
    random_values=np.random.uniform(size=the_sample)
    result=[]
    for i in random_values:
        if i >= cumulative[1].min():
            res = minimize_scalar(spline_fit_function, bounds=[cumulative[0].min(),np.log(average+1),cumulative[0].max()], args=(i,inv_cumulative,0), method='Golden')
        #result.append(np.exp(res.x)-1)
            result.append(res.x)
        else:
            result.append(0.0)

    return np.array(result)#}}}



if __name__ == '__main__':
	params = [i[1:5] for i in readFile(av[1])]
	res = []
	dump = open(av[2],'w')
	errC = 0
#	print mrna_sample_probability(l=1.0,g=1.0,mu=8.0,d=0.1, the_sample=1000)
	for i in params:
		k0,k1,v0,d0 = i[0],i[1],i[2],i[3]		
		draw = mrna_sample_probability(l=k0,g=k1,mu=v0,d=d0, the_sample=1000)
		try:
			out = [np.mean(draw),(np.std(draw)/np.mean(draw))**2,mrna_average(k0=k0,k1=k1,v0=v0,d0=d0),mrna_variance(k0=k0,k1=k1,v0=v0,d0=d0)/(mrna_average(k0=k0,k1=k1,v0=v0,d0=d0))**2]
			dump.write('\t'.join([str(k) for k in out])+'\n')
			dump.flush()
		except: 
			errC += 1
			print 'You had bad luck %i times today. What a pity!' %(errC)
