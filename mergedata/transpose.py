import sys

def readData(inname):
	data = []
	f = open(inname,'r')
	for l in f:
		parts = l.strip().split('\t')
		data.append(parts)
	f.close()
	return data

def trans(data):
	tdata = [[data[i][j] for i in range(len(data))] for j in range(len(data[0]))]
	return tdata

def writeData(outname,data):
	f = open(outname,'w')
	for vec in data:
		f.write('%s\n' % '\t'.join(vec))
	f.close()

if __name__ == '__main__':
	data = readData(sys.argv[1])
	tdata = trans(data)
	writeData(sys.argv[2],tdata)
