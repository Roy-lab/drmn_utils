import sys
import argparse

parser = argparse.ArgumentParser('Make per cluster regression matrix')
parser.add_argument("--config",    type=str, help="DRMN config file",            required=True)
parser.add_argument("--order",     type=str, help="Cell line orders",            required=True)
parser.add_argument("--clusters",  type=str, help="transitioning gene clusters", required=True)
parser.add_argument("--outprefix", type=str, help="Prefix of output",            required=True)
args = parser.parse_args()

def readOrder(inname):
	times = []
	f = open(inname,'r')
	for l in f:
		times.append(l.strip())
	f.close()
	return times

def readCluster(inname):
	clusters = {}
	f = open(inname,'r')
	for l in f:
		parts = l.strip().split('\t')
		g = parts[0]
		c = parts[1]
		genes = clusters.get(c,set([]))
		genes.add(g)
		clusters[c] = genes
	f.close()
	return clusters

def readExp(inname,suf):
	cnt = 0
	edata = {}
	f = open(inname,'r')
	for l in f:
		cnt += 1
		if cnt == 1:
			continue
		parts = l.strip().split('\t')
		g = parts[0].replace(suf,'')
		edata[g] = float(parts[1])
	f.close()
	return edata

def readFeat(inname,suf):
	cnt = 0
	fdata = {}
	tfs = set([])
	f = open(inname,'r')
	for l in f:
		cnt += 1
		if cnt == 1:
			continue
		parts = l.strip().split('\t')
		tf = parts[0]
		tg = parts[1].replace(suf,'')
		tfs.add(tf)
		fdata[(tf,tg)] = float(parts[2])
	f.close()
	return (tfs,fdata)

def writeData(outname,times,clusters,allexp,allfeat,alltfs):
	tfs = list(alltfs)
	tfs.sort()
	#times = allexp.keys()
	#times.sort()
	for c in clusters:
		genes = list(clusters[c])
		genes.sort()
		print('Writing for cluster%s...' % c)
		f = open('%s_%s.txt' % (outname,c),'w')
		f.write('Gene\tExp\t%s\n' % '\t'.join(tfs))
		for t in times:
			edata = allexp[t]
			fdata = allfeat[t]
			for g in genes:
				gname = '%s_%s' % (g,t)
				vec = ['%f' % fdata.get((tf,g),0) for tf in tfs]
				f.write('%s\t%f\t%s\n' % (gname,edata.get(g,0),'\t'.join(vec)))
		f.close()

if __name__ == '__main__':
	times = readOrder(args.order)
	allexp = {}
	allfeat = {}
	alltfs = set([])
	f = open(args.config,'r')
	cnt = 0
	for l in f:
		cnt += 1
		if cnt == 1:
			prefix = l.strip()
			continue
		parts = l.strip().split('\t')
		t = parts[0]
		suf = '_%s' % t
		print('reading expression for %s...' % t)
		edata = readExp(prefix+'/'+parts[2],suf)
		print('reading features for %s...' % t)
		(tfs,fdata) = readFeat(prefix+'/'+parts[3],suf)
		allexp[t] = edata
		allfeat[t] = fdata
		alltfs = alltfs.union(tfs)
	f.close()
	print('reading clusters...')
	clusters = readCluster(args.clusters)
	writeData(args.outprefix,times,clusters,allexp,allfeat,alltfs)
