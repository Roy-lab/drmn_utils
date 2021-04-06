import sys

def readMot(inname):
	cmap = {}
	f = open(inname,'r')
	for l in f:
		#chr1	3004298	3004313	+	8.830266
		parts = l.strip().split('\t')
		c  = parts[0]
		p1 = int(float(parts[1]))
		p2 = int(float(parts[2]))
		ID=parts[3];
		s  = float(parts[4])
		ss = parts[5]
		t  = cmap.get(c,[])
		t.append((p1,p2,s,ss,ID))
		cmap[c] = t
	f.close()
	return cmap

def readGenes(inname):
	gmap = {}
	f = open(inname,'r')
	for l in f:
		#ENSMUST00000166088	chr10	75032585	75032585
		parts = l.strip().split('\t')
		g  = parts[0]
		c  = parts[1]
		p1 = int(float(parts[2]))
		p2 = int(float(parts[3]))
		strand=parts[4]
		gs = gmap.get(c,[])
		gs.append((g,p1,p2,strand))
		gmap[c] = gs
	f.close()
	return gmap

# take second element for sort
def takeFirst(elem):
    return elem[0]

def takeWindow(elem,w_up,w_down):
	if elem[3] == '+':
    		return elem[1]-w_up
	else:
		return elem[2]-w_down

def writeTSS(outname,cmap,gmap,w_up,w_down):
	f = open(outname,'w')
	genecnt = {}
	for (c) in gmap:
		if (c) not in cmap:
                        continue
		#sort genes by starting edge of mapping window for this gene
		gs = gmap[c]
		gs.sort(key=lambda x:takeWindow(x,w_up,w_down))
		#sort by starting coordinate of motif.
		cs = cmap[c]
		cs.sort(key=takeFirst)
		start=0
		#loop over motif entries
		for (mp1,mp2,v,ss,ID) in cs:
			#loop over gene entries
			for gI in range(start,len(gs)):
				print start
				#get gene information 
				(g,p1,p2,gstrand)=gs[gI]
				#if + strand gene check
				if gstrand == '+':
					if p1+w_down < mp1:
						start=gI
						continue
						#start=gI
					if p1-w_up > mp2:
						break
					if mp1 >= p1-w_up and mp2 <= p1+w_down:
						f.write('%s\tCONVERT\tTSS_%s_TSS%d\t%d\t%d\t%f\t%s\t.\t%s;%s\n' % (c,g,1,mp1,mp2,v,"+",ID,g));
				#if - strand gene check in slightly differnet way, determining the window coordinates differently 
				else:
					if p2+w_up < mp1:
						start=gI
						continue
						#start=gI
					if p2-w_down > mp2:
						break
					if mp1 >= p2-w_down and mp2 <= p2+w_up:
						f.write('%s\tCONVERT\tTSS_%s_TSS%d\t%d\t%d\t%f\t%s\t.\t%s;%s\n' % (c,g,1,mp1,mp2,v,"+",ID,g));
	f.close()

if __name__ == '__main__':
	cmap = readMot(sys.argv[1])
	gmap = readGenes(sys.argv[2])
	w_up = int(sys.argv[3])
	w_down = int(sys.argv[4])
	net = writeTSS(sys.argv[5],cmap,gmap,w_up,w_down)
