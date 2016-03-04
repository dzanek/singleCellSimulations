import stochpy as sp
import random as r
import numpy as np
import sys


def readParams(fname):
	data = open(fname,'r').readlines()
	data = [[float(i) for i in k.split()] for k in data if k[0] != '#']
	return data


def randGenes(genes,gnumber):
	data = r.sample(genes,gnumber)
	return data

def doStochSim(paramList,steps,smodel,ncell): 
	#here it takes: list of 4 params, number of time units for sim.length., name of model (obsolate), number of cells to sim.
	
	
        datList = []
	print ncell			
	for i in range(int(ncell)):				
		sim.ChangeInitialSpeciesAmount('mRNA',int(paramList[5]))
	        sim.ChangeParameter("kdeg",paramList[4])	
	        sim.ChangeParameter("kon",paramList[1])
	        sim.ChangeParameter("koff",paramList[2])	
		sim.ChangeParameter("ksyn",paramList[3])
	        sim.DoStochSim(mode='Time',end=steps)
	        sim.GetTrajectoryData()
		valAtTime = sim.data_stochsim.getAllSimData()[-1][-1] 
		datList.append(valAtTime)
		print '>>><<<\n',valAtTime,'\t',valAtTime-paramList[5],'\t',paramList[5],'\n>>><<<\n'
	mean = np.mean(datList)
	std = np.std(datList)
	'''
	wTimeAvg0 = np.mean([i[0] for i in wTime])
	wTimeAvg1 = np.mean([i[1] for i in wTime])
	wTimeAvg2 = np.mean([i[2] for i in wTime])
	wTimeAvg3 = np.mean([i[3] for i in wTime])
	'''
	print 'Simulated: ',mean,'\t',(std/mean)**2
	sqvT = (paramList[5]*(1+((paramList[2]*paramList[3])/(((paramList[1]+paramList[2])**2)+(paramList[4]*(paramList[1]+paramList[2]))))))/(paramList[5]**2)
	print 'Analytical: ',paramList[5],'\t',sqvT
        tr = [str(mean),str((std/mean)**2),paramList[5],sqvT,str('\t'.join([str(i) for i in paramList[1:5]]))]+datList
	return tr #here is: mean, CVsq for sim and calc (4 colls), parameters (4 colls), single cell expression levels 


if __name__ == '__main__':
	params = readParams(sys.argv[1])	#read param file from genGenes.py
	sample = params[int(sys.argv[-1]):int(sys.argv[-1])+2000]# it's good to take slice of parameter// also You can sample it randGenes(params,int(sys.argv[2]))
	results = []
	dump = open(sys.argv[6]+sys.argv[-1],'w') #name for output with slice marker	
	counter = 0			
        sim = sp.SSA(Method='TauLeaping')	
        sim.Model('Burstmodel.psc')
	for i in sample:
		i = ['nan'] + i	#to add first column was simplier than change rest of code ;)
		counter += 1
		print counter
		dump.write('\t'.join([str(v) for v in doStochSim(i,float(sys.argv[3]),sys.argv[5],sys.argv[4])])+'\n') #here You run simulation function
		dump.flush()	#i like having results in file when they are ready

