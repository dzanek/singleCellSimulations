import rf
import sys
import numpy as np
import matplotlib.pylab as pl
import wf
import scipy.special as ssp
import scipy.stats as ss
import random as r

def pcr_binomial(start_n,p, nr_pcr=15):
	#takes input molecules, probability and cycles.
    result=[]
    for i in range(1):
        final_n=start_n
        for i in range(nr_pcr):
            final_n+= np.random.binomial(final_n,p)
        result.append(final_n)
    return np.array(result)[0]


def calcNoise(m0,p,n):

	mean = m0*(1+p)**n
	std = m0*((1-p)/(1+p))*((1+p)**(2*n)  - (1+p)**n )
	return (mean,std)

def noneg(q):
	if q>1: return int(q)
	else: return 0

	'''
def fragmentation(length):
		
	inLen = length
	fragNo = 1
	while length>10:
		length = length-ss.weibull_min.rvs(np.log10(inLen),scale=200)
		fragNo += 1
	return fragNo


	'''	
def pcr_sample_binomial(start_n,p, nr_pcr, the_sample=1):
    result=[]
    for i in range(the_sample):
        final_n=start_n
        for i in range(nr_pcr):
            final_n+= np.random.binomial(final_n,p)
        result.append(final_n)
    return np.array(result)


def mkLib(inp,p,nr_pcr):
	
	if inp == 0: return 0	#no expression no seq
	if inp<1: inp = np.random.binomial(1,inp*0.1)	#if expression id <1 just make probability lower
	else: inp = np.random.binomial(inp,0.3) #more abudant fragments more likely comes to pcr
	if inp == 0: return 0			#again - if no fragment came to PCR effeciency is 0
	else: return pcr_binomial(inp,p,nr_pcr)
	'''
	#fragmentation may need including. 
	for copy in range(int(inp)):
		fromCopy = 0
		fragNo = fragmentation(length)
		for frag in range(int(fragNo)):
			#if r.uniform(0,1) < 0.1: fromCopy += pcr_binomial(1,p,nr_pcr)			
			fromCopy = pcr_binomial(1,p,nr_pcr)
		total += fromCopy
	'''
	return total






expr = [i[int(sys.argv[2])] for i in rf.readFile(sys.argv[1])[1:] if i[int(sys.argv[2])] != 0]
#reads  file with expression values from certain column (files I used were various)

for p in np.linspace(0.05,0.25,10): #probe few PCR effeciency levels
	p=round(p,2)
	results = []
	pC = 34
	for i in expr:	#for each gene
		fragno = 1 #obsolate modification accounting for fragmentation //ss.lognorm.rvs(0.7,loc=0,scale=2600)
		results.append([v for v in [i,fragno*mkLib(i,p,pC),fragno*mkLib(i,p,pC),fragno*mkLib(i,p,pC),fragno*mkLib(i,p,pC),fragno*mkLib(i,p,pC),fragno*mkLib(i,p,pC)]])
		#appends with number of: expressed, techrepl1, techrepl2, more techrepls 
	tocor = [i for i in results if 0.0 not in i[1:3]]	#gets rid of zeroes for log10
	print p,len(results), len([i for i in results if i[2]!=0 and i[1]!=0 ]), len([i for i in results if i[2]==1 and i[1]==1 ]), np.corrcoef([np.log10(i[1]) for i in tocor],[np.log10(i[2]) for i in tocor])[0][1] 
	cc =  np.corrcoef([np.log10(i[1]) for i in tocor],[np.log10(i[2]) for i in tocor])[0][1]
	alfa = beta = 0
	#part for ploting two random replicates
	while alfa == beta:
		alfa, beta = r.randint(1,6), r.randint(1,6)
	if alfa != beta: pl.plot([i[alfa] for i in results],[i[beta] for i in results],'o',label='%s %i %i' %(p,alfa,beta))
	pl.xscale('log')
	pl.yscale('log')
	pl.xlabel('Corr.coef. ='+str(cc))
	pl.savefig('pcrNoise_flux%s.png' %(p))
	pl.close()
	dump = open ('%s.dat' %(p),'w')
	dump.write('ex \t a \t b \t c \t d \t e \t f \t li')
	dump.write('\n')
	for i in results:
		dump.write('\t'.join([str(k) for k in i])+'\n')
	dump.close()






