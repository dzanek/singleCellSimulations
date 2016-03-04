def writeFile(fname,data):
	dump = open(fname,'w')
	for k in data:
		dump.write('\t'.join([str(i) for i in k])+'\n')
		dump.flush()

