#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Framework.H"

/*
* Given per-base signal (from genomecov) and a set of regions with TSS names, 
* motif names, and gene names.
* Each gene (last column) can have multiple TSSes (third column), which in turn
* may have multiple regions (first column).
* We sum up all regions with the same TSS name.
* For each gene, we take the max TSS.
*/


Framework::Framework(bool myPrintProfile, int myWindowSize)
{
		total=0;
		minVal=0;
		dist=0;
		motifNet=true; // always print motif network
		printProfile=myPrintProfile;
		windowSize=myWindowSize;
		doNormalize=false;
		normF=0;
}

Framework::~Framework()
{
}

/*
* Reads the gtf-lite format. 
*
* Example:
* chr3    CONVERT TSS_GENE1_TSS1    99535755        99537755        .       +       .       Motif1;ENST00000496116
* chr20   CONVERT TSS_GENE2_TSS1    57225368        57227368        .       +       .       Motif2;ENST00000496117
* chr1    CONVERT TSS_ENSG00000223972_ENST00000456328.2   11226   11240   518     +       .       JASPAR_300;ENSG00000223972
*
* (0) Chromosome
* (1) Not used
* (2) TSS name. We will sum up all motif regions for this TSS. Looks like: TSS_ENSG00000227232_ENST00000537342.1
* (3) Start position (if + strand), end position (if - strand)
* (4) End position (if + strand), start position (if - strand)
* (5) Not used
* (6) Strand
* (7) Not used
* (8) Motif;Gene name. We will end up reporting the max value over all TSSes associated with this motif AND gene.
*/

int 
Framework::readTSS(const char* aFName)
{
	int success=0;
	ifstream inFile(aFName);
	char buffer[4096];
	cout << "Window size" << windowSize << endl;

	while(inFile.good())
	{
		inFile.getline(buffer,4095);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string chromosome;
		string tssgene; // make from 8 and 2. looks like JASPAR_300_TSS_ENSG00000223972_ENST00000456328.2, I'm not sure why this exists??
		string tssname; // field 2 -- name for this TSS and gene
		string motifGenePair; // field 8 (actually a pair)
		char strand;
		int position_b;
		int position_e;
		int position;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				chromosome.append(tok);
			}
			else if (tokCnt==2)
			{
				tssname.append(tok);
			}
			else if(tokCnt==3)
			{
				position_b=atoi(tok);
				//flipping
				//position_e=atoi(tok);
			}
			else if(tokCnt==4)
			{
				position_e=atoi(tok);
				//flipping
				//position_b=atoi(tok);
			}
			else if(tokCnt==6)
			{
				strand=tok[0];
			}
			else if(tokCnt==8)
			{
				char* pos=strchr(tok,'.');
				if(pos!=NULL)
				{	
				//	*pos='\0';
				}
				pos=strchr(tok,' ');
				if(pos!=NULL)
				{
					tok=pos+1;
				}
				pos=strchr(tok,' ');
				if(pos!=NULL)
				{	
					*pos='\0';
				}
				motifGenePair.append(tok);
				// stick the motifGenePair on the tss also. I think there was a reason for this, but IDK anymore.
				tssgene.append(tok);
				tssgene.append("_");
				tssgene.append(tssname);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		if(motifGenePair.length()<=0)
		{
			continue;
		}
		if(chromosome.find("MT")!=std::string::npos)
                {
			continue;
		}
		map<string,Framework::MotifTSS*>* edgeSet=NULL;
		map<int,Framework::MotifTSS*>* tsses=NULL;
		
		// for this chromosome, do we have a set of potential edges yet?
		if(tssSet.find(chromosome)==tssSet.end()) 
		{
			edgeSet=new map<string,Framework::MotifTSS*>;
			tssSet[chromosome]=edgeSet;
			tsses=new map<int,Framework::MotifTSS*>;
			tsscoordSet[chromosome]=tsses;
		}
		else
		{
			edgeSet=tssSet[chromosome];
			tsses=tsscoordSet[chromosome];
		}
		if(strand=='+')
		{
			position=position_b;
		}
		else	
		{
			position=position_e;
		}
		
		// make/get gene struct: store tssgene
		Framework::MotifTSS* tss=NULL;
		if(edgeSet->find(tssgene)==edgeSet->end())
		{
			tss=new Framework::MotifTSS;
			(*edgeSet)[tssgene]=tss;
			tss->edgeTriplet.append(tssgene); // store tss/motif here
			tss->strand=strand; // add strand info
			tss->motifGenePair.append(motifGenePair); // add gene/motif here
		}
		else
		{
			tss=(*edgeSet)[tssgene];
		}
		
		// sanity check -- strand the same as the one we already stored?
		if (strand != tss->strand)
		{
			cerr << "How weird -- found both sense and antisense for same motif/TSS?? " << motifGenePair << endl;
		}
		
		// fill or append gene TSS bin information
		tss->tsses[position]=0; // store this instance for this TSS
		(*tsses)[position]=tss; // and store this position to correspond to this TSS
		if(strand=='+')
		{
			tss->beginnings.push_back(position_b);
			tss->endings.push_back(position_e);
		}
		else // swap beginning/end b/c antisense 
		{
			tss->beginnings.push_back(position_e);
			tss->endings.push_back(position_b);
		}	
		// make sure this motif instance fits within our defined window size
		double regSize=abs(position_b-position_e);
		if (printProfile && (regSize > 2*windowSize))
		{
			cout << "Sorry, motif instance of size " << regSize
			<< " is too long to fit into a window [-" << windowSize << ", " << windowSize << "]" << endl;
			success=1;
			break; 
 		}
	}
	inFile.close();
	return success;
}

/*
* Reads chromosome sizes info and populates a massive data structure with room 
* to accept the read info.
*/ 
int
Framework::readChromosomeSizes(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string chromosome;
		int size;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				chromosome.append(tok);
			}
			else if(tokCnt==1)	
			{
				size=atoi(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		if(chromosome.find("MT")==std::string::npos)
		{
			dist+=(double)size;
			chromosomeSizes[chromosome]=size;
		}
	}
	inFile.close();
	return 0;
}
	


int 
Framework::readSignals(const char* aFName)
{
	readData(aFName,signal,signalAvg);
	return 0;
}


/**
* Normalize stored signal by total and multiply by provided scalar.
*/
int 
Framework::normalize(double scalar)
{
	double normF=scalar*total/dist;
	cout << "normalization by " << normF*scalar << endl;
	/*cout << "from " << total << " over " << dist << endl;
	for(map<string,vector<double>*>::iterator cIter=signal.begin();cIter!=signal.end();cIter++)
	{
		string chromosome=cIter->first;
		vector<double>* values_fg=signal[cIter->first];
		for(int i=0;i<values_fg->size();i++)
		{
			double normVal=log(((*values_fg)[i]/normF*scalar)+1);
			(*values_fg)[i]=normVal;
		}
	}*/
	doNormalize=true;
	return 0;
}

int 
Framework::readData(const char* fName,map<string,vector<double>* >& data, double& globalAvg)
{
	ifstream inFile(fName);
	char buffer[1024];
	int cnt=0;
	double totalBins=0;
	globalAvg=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer," \t");
		int tokCnt=0;
		string chr;
		int begin;
		int end;
		double val;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				chr.append(tok);
			}
			else if(tokCnt==1)
			{
				begin=atoi(tok);
			}
			else if(tokCnt==2)
			{
				end=atoi(tok);
			}
			else if(tokCnt==3)
			{
				val=atof(tok);
			}
			tok=strtok(NULL," \t");
			tokCnt++;
		}
		// add pseudocount to avoid 0? not now.
		// seems like it isn't doing what we want it to do.
		//val=val+1e-3; // TODO I don't think this is doing what we want it to do
		
		if(chr.find("MT")!=std::string::npos)
		{
			continue;
		}

		total+=(end-begin)*val;

		vector<double>* values=NULL; // allocate entire chromosome
		
		if(data.find(chr)==data.end())
		{
			if(chromosomeSizes.find(chr)==chromosomeSizes.end())	
			{
				cout <<"Chromosome sizes are not read in or chromosome " << chr << " does not exist" <<endl;
			}
			int chrSize=chromosomeSizes[chr];
			values=new vector<double>(chrSize);
			for(int i=0;i<chrSize;i++)
			{
				(*values)[i]=0;
			}
			data[chr]=values;
			cout <<"Created chromosome " << chr<< endl;	
		}
		else
		{
			values=data[chr];
		}
		//Now store the value for every base pair in the count region
		for(int i=begin;i<end;i++)
		{	
			(*values)[i]=val;
		}
		cnt++;
		if(cnt%500000==0)
		{
			cout <<".";
			cout.flush();
		}
		globalAvg=globalAvg+val;
	}	
	inFile.close();
	totalSignal=globalAvg; // store total
	globalAvg=globalAvg/cnt;
	cout << endl <<"Read "<< cnt << " values " << " Avg: " << globalAvg  << endl;
	cout << endl <<"Global total was " << totalSignal << endl;
	return 0;
}


/*
* Gets the SUM signal value for each region associated with a name.
* Prints to file with name aFName.
* For fun we also print tne number of instances in an additional column.
* If motifnet=true, then we assume name is format motif;gene and convert ; to tab.
* We ALSO print out the per-
*/
int 
Framework::getSignalInRegion(const char* aFName) //, bool motifnet, bool printProfile)
{
	// make signal file
	char sigFName[1024];
	sprintf(sigFName,"%s_signalProfile.txt",aFName);
	ofstream sigFile; 

	// only open file if we are printing the profile
	if (printProfile)
	{
		sigFile.open(sigFName);
	}

	char avgFName[1024];
	sprintf(avgFName,"%s_aggregate.txt",aFName);
	ofstream avgFile(avgFName);
	
	//Go over each chromosome and output the avg value for a gene.
	// Take the SUM over possible values
	cout << "Now going over " << tssSet.size() << " chromosomes" << endl;
	for(map<string,map<string,Framework::MotifTSS*>*>::iterator cIter=tssSet.begin();cIter!=tssSet.end();cIter++)
	{
		map<string,Framework::MotifTSS*>* edgeSet=cIter->second;
		cout << "Working on chromosome " << cIter->first << " with " << edgeSet->size() << " motif-TSS combinations" << endl;
		if(signal.find(cIter->first)==signal.end())
		{
			cout << "Skipping chromosome " << cIter->first << endl;
			continue;
		}

		int success=processChromosome(cIter->first, avgFile, sigFile);
		if (success != 0)
		{
			cerr << "Encountered a problem on chromosome " << cIter->first << endl;
		}
		// clear this chromosome's signal
		signal[cIter->first]->clear(); // clears up input data for this chromosome
		delete signal[cIter->first];
	}
	// clear out the signals

	//close the open files
	if (printProfile)
	{
		sigFile.close();
	}

	avgFile.close();
	return 0;
} // end get signal
	

/*
* Aggregates signal for one chromosome and writes to open files.
* Steps:
*  
*/
int
Framework::processChromosome(string chromosome, ofstream& avgFile, ofstream& sigFile)
{	
	map<string,Framework::MotifTSS*>* edgeSet=tssSet[chromosome];

	vector<double>* svals=signal[chromosome]; // stores values for the whole chromosome

	// store for this chromosome
	// { gene/motif : {gene/motif/tss : aggsignal}}
	map<string, map<string, double>* > storedSignal;

	// hold onto number of motif instances in a tss bin
	// { gene/motif : { gene/motif/tss : numMotifInstances }}
	//map<string, map<string, int>* > storedNumRegions;
	// key: gene/motif/tss; val: number of motif instancesin TSS
	map<string, int> storedNumRegions;
	

	// the first thing we have to do is get scores and profiles for each motif/tss/gene edge
	// key: edge name (motif-gene-tss triple) ; val: the edge object(one line in the input file)
	// this populates the storedSignal and storedNumRegions data
	int popOK=populateTSSRecords(chromosome, storedSignal, storedNumRegions);


	// The second thing we do is loop over the motif-gene edges to figure out which is the best TSS
	// and print that info out to our file(s).
	int reportOK=reportTSSRecords(chromosome, storedSignal, storedNumRegions, avgFile, sigFile);

	// clean up! 
	// storedNumRegions doesn't need to be cleaned up - in score
	// need to clean: storedSignal
	clearStoredSignal(storedSignal);

	return 0;
} 

/**
* collect scores and per-bp profiles for each motif/gene/tss edge
* Populates: storedSignal, storedNumRegions
*/
int
Framework::populateTSSRecords(
	string chromosome, 
	map<string, map<string, double>* >& storedSignal, 
	map<string, int>& storedNumRegions)	
{
	// grab the TSS records to populate
	map<string,Framework::MotifTSS*>* edgeSet=tssSet[chromosome];

	// grab the signal data for this chromosome
	vector<double>* svals=signal[chromosome]; // stores values for the whole chromosome

	for(map<string,Framework::MotifTSS*>::iterator gIter=edgeSet->begin();gIter!=edgeSet->end();gIter++)
	{
		Framework::MotifTSS* g=gIter->second;

		double sumForMotifTSS=0; // saves the sum value for this motif/gene/TSS instance
		int numRegions=0;

		//loop over all motif instances within this particular TSS (beginnings) 
		for(int i=0; i<g->beginnings.size(); i++)
		{
			int binstart=g->beginnings[i];
			int binend=g->endings[i];
			// sanity check against chrom sizes -- if bin goes off the end, make it stop at the last position
			binend=min(binend, chromosomeSizes[chromosome]);
			
			double mean=0;
			//vector<double> vset;
			for(int b=binstart;b<=binend;b++)
			{
				int coord_key=b;
				mean=mean+(*svals)[coord_key];	
			}

			if (mean==0)
			{
				continue;
			}
			

			// take avg counts over the motif instance.
			mean=mean/(binend-binstart+1);
			sumForMotifTSS+=mean; // add this motif instance
			numRegions++;
			//cout << "\t" << g->edgeTriplet << "\t" << mean << endl;

		} // end loop over beginnings

		if(sumForMotifTSS<=minVal)
		{
			cout <<"This region" << g->edgeTriplet << " really has no signal." << endl;
			continue;
		}
		
		// Store sum AND each instance's per-bp signal for this motif/tss/gene comb
		// motifGenePair is the motif/gene combo, from which we need to pick the best TSS
		// edgeTriplet is the motif/tss/gene combo
		if (storedSignal.find(g->motifGenePair) == storedSignal.end())
		{
			storedSignal[g->motifGenePair]=new map<string,double>();
		}
		// for this motif/gene/TSS, store sum and number of instances
		(*storedSignal[g->motifGenePair])[g->edgeTriplet]=sumForMotifTSS;
		storedNumRegions[g->edgeTriplet]=numRegions;
	} // end loop over triplets

	return 0;
}

/**
* Goes through the populated records. Chooses best TSS for each motif-gene pair and reports 
* the aggregated signal. 
* If desired (global var printProfile set to true), prints out the per-bp signals for each motif instance.
* Prerequisite: populateTSSRecords has been run
*/
int
Framework::reportTSSRecords(string chromosome, 
	map<string, map<string, double>* >& storedSignal, 
	map<string, int>& storedNumRegions,
	ofstream& avgFile,
	ofstream& sigFile)
{
	// for each motif-gene pair, find best TSS and report results.
	for (map<string, map<string,double>*>::iterator gIter=storedSignal.begin(); gIter!=storedSignal.end(); gIter++)
	{
		string motifGenePair=gIter->first;
		string bestTSS;
		double maxSignal=-1;
		map<string,double>* tssMap=gIter->second;
		// loop over possible TSSes for this gene; find the best.
		for (map<string,double>::iterator tssIter=tssMap->begin(); tssIter!=tssMap->end(); tssIter++)
		{
			if (tssIter->second > maxSignal)
			{
				maxSignal=tssIter->second;
				bestTSS=tssIter->first;
				//cout << "Updating " << motifGenePair << " " << bestTSS << " " << maxSignal << endl;
			}
		}


		// Now print the aggregate motif-gene signal to file. if motifNet is true, we'll expand motif;gene into motif\tgene.
		// edgeName is the motif/gene combo...
		string edgeName=motifGenePair; 
		if (motifNet)
		{
			edgeName=motifGenePair.replace(motifGenePair.find(";"),1,"\t");
		}		

		if (maxSignal > 0)
		{
			// sum avg signal for a gene, number of instances summed up, best TSS 
			int numRegions=storedNumRegions[bestTSS];
			//SK:update value to the log ratio relative to the global mean if normalization is requested with -n 
			if(doNormalize)
			{
				double normVal=log((maxSignal/normF)+1);
				maxSignal=normVal;	
			}
			avgFile << edgeName <<"\t" << maxSignal << "\t" << numRegions << "\t" << bestTSS << endl;			
			// print the profiles if desired
			if (printProfile)
			{
				printInstanceProfiles(sigFile, chromosome, bestTSS);
			}
		}
	} // end loop over gene-motif 

	return 0;
}

/*
* Clears the memory from the entries in storedSignal.
* Each key is one motif->gene edge; the values are the possible TSSes that could
* be used to score the edge.
*/
int
Framework::clearStoredSignal(map<string, map<string, double>* >& storedSignal)
{
	// for each motif-gene pair
	for (map<string, map<string,double>*>::iterator gIter=storedSignal.begin(); gIter!=storedSignal.end(); gIter++)
	{
		string motifGenePair=gIter->first;
		map<string,double>* tssMap=gIter->second;
		tssMap->clear();
	} // end loop over gene-motif 
	storedSignal.clear(); //  fully clear out the keys/values
	return 0;
}


int 
Framework::printInstanceProfiles(ofstream& sigFile, string chromosome, string bestTSS)
{
	//cout << "Now printing per-BP signal for all motif instances in " << bestTSS << endl;

	//loop over all motif instances within this particular TSS (beginnings) 
	Framework::MotifTSS* triple=(*tssSet[chromosome])[bestTSS];

	// access the values for this chromosome
	vector<double>* svals=signal[chromosome];

	for(int i=0; i<triple->beginnings.size(); i++)
	{

		// get expanded window around center
		int instanceStart=triple->beginnings[i];
		int instanceEnd=triple->endings[i];
		int instanceCenter=instanceStart + (instanceEnd-instanceStart)/2;

		// name this instance for the actual coordinates
		char buffer[1024]; // chr_start_end|motif_gene_tss
		sprintf(buffer, "%s_%d_%d|%s", chromosome.c_str(), instanceStart, instanceEnd, bestTSS.c_str());
		string rowName(buffer);

		int binstart=instanceCenter-windowSize;
		int binend=instanceCenter+windowSize+1;

		sigFile << rowName;		
		for(int b=binstart;b<=binend;b++)
		{
			// sanity check against chromosome size: print NaN if we go off the end.
			if (b < 0 || b >= chromosomeSizes[chromosome])
			{
				sigFile << "\tNaN";
			}
			else if(!doNormalize)
			{
				sigFile << "\t" << (*svals)[b]; // print stored signal value
			}
			else if(doNormalize)
			{
				double normVal=log(((*svals)[b]/normF)+1);
				sigFile << "\t" << normVal;
			}
		}
		sigFile << endl; // end this row
	} // end loop over beginnings

	return 0;
}


// updated by DC for motifnet, and added profile printing back
// added seq depth norm option
int
main(int argc, const char** argv)
{
	string USAGE="Usage: aggregateSignalToSum tss chromosomesizes signal outfilePrefix [-p/--printProfile <windowSize>, -n/--normalize <scalar>]\n"
		"If -p/--printProfile option, will print per-base-pair signal for motif instances. Use with caution as it will create huge files.\n"
		"\tThis requires a parameter for defining a uniformly sized window around each motif center: center +/- window.\n"
		"\tYou need to set the window size at least as wide as half of your largest motif.\n"
		"If -n/--normalize <scalar>, will normalize to a sequencing depth of 1 million reads.";

	if(argc<5 || argc > 8)
	{
		cout << USAGE << endl;
		return 2;
	}

	bool motifNet=true;
	bool printProfile=false;
	bool doNormalize=false; 
	int windowSize=0;
	double normScalar=0; 

	const char* tssFile; //arg 1
	const char* sizeFile; // arg 2
	const char* signalFile; // arg 3
	const char* outputFile; // arg 4

	int a=1; // argument counter (we'll skip over options)
	// i counts over things in argv
	for (int i=1; i<argc; i++)
	{
		string argstring(argv[i]);
		bool isValid=true; // did we match an option or argument?

		// either long option (--blah) or short option alone (-b); value is in next token.
		if (argstring.find("--")==0 || (argstring.find("-")==0 && argstring.size()==2))
		{
			cout << "Got option, value in next token: " << argstring << " " << argv[i+1] << endl;
			if (argstring.find("-p") != string::npos)
			{
				windowSize=atoi(argv[i+1]);
				cout << "Got window size " << windowSize << " from next token " << argv[i+1] << endl;
				printProfile=true;
				i++;
			}			
			else if (argstring.find("-n") != string::npos)
			{
				normScalar=atof(argv[i+1]);
				cout << "Got normalization factor " << normScalar << " from next token " << argv[i+1] << endl;
				doNormalize=true;
				i++;
			} 
			else
			{
				isValid=false;
			}
		} 
		// otherwise: we're looking at the short option (-p/-n) with the value in the same token.
		else if (argstring.find("-")==0 && argstring.size()>2)
		{
			if (argstring.find("-p")==0)
			{
				// atoi returns 0 if integer not at beginning of string
				windowSize=atoi(argstring.substr(2).c_str()); 
				printProfile=true;
			} 
			else if (argstring.find("-n")==0)
			{
				// atof returns 0 if integer not at beginning of string
				normScalar=atof(argstring.substr(2).c_str()); 
				doNormalize=true;
			}	
			else
			{
				isValid=false;
			}
		}
		// argument parsing
		else
		{
			isValid=true; // assume valid unless not
			switch(a)
			{
				case 1: tssFile=argv[a]; break;
				case 2: sizeFile=argv[a]; break;
				case 3: signalFile=argv[a]; break;
				case 4: outputFile=argv[a]; break;
				default: isValid=false; break;
			}
			// proceed to next argument
			a++;
		}
		
		if (!isValid)
		{
			cout << "Illegal argument: " << argv[i] << endl;
			cout << USAGE << endl;
			return(2);
		}
	}

	// double-check settings for options
	if (printProfile && windowSize <=0)
	{
		cout << "Window size for profile printing must be > 0. You provided: " << windowSize << endl;
		return 2;
	} 
	else if (doNormalize && normScalar <=0) 
	{
		cout << "Normalization factor must be > 0. You provided: " << normScalar << endl;
		return 2;
	}

	Framework fw(printProfile, windowSize);
	int tssOK=fw.readTSS(tssFile);
	if (tssOK != 0)
	{
		cerr << "Problem reading motif instance file." << endl;
		return 1;
	}
	fw.readChromosomeSizes(sizeFile);
	fw.readSignals(signalFile);
	if (doNormalize)
	{
		fw.normalize(normScalar);
	}
	fw.getSignalInRegion(outputFile); 
	return 0;
}
	
