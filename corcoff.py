import numpy as np
import sys


plik = open(sys.argv[1],'r').readlines()
plik = [i.split() for i in plik]
plik = [[float(k) for k in i[:8]] for i in plik if np.isnan(float(i[1]))==np.isnan(float(i[3]))==False]
if len(sys.argv) > 2 and sys.argv[2] == 'both':
	calc = [i[0] for i in plik]
	ther = [i[1] for i in plik]
	a = np.corrcoef(calc,ther)
	print a[0][1]

if len(sys.argv) > 2 and sys.argv[2] == 'log':
	calc = [np.log10(i[1]) for i in plik]
	ther = [np.log10(i[3]) for i in plik]
	a = np.corrcoef(calc,ther)
	print sys.argv[1],sys.argv[1][4:6],a[0][1],'\tlog'
	
calc = [i[1] for i in plik]
ther = [i[3] for i in plik]
a = np.corrcoef(calc,ther)
print sys.argv[1],sys.argv[1][4:6],a[0][1]
