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
