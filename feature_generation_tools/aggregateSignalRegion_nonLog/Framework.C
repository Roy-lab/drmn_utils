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
Framework::readRegions(const char* aFName)
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
		string name;
		char strand='+';
		int position_b;
		int position_e;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				chromosome.append(tok);
			}
			else if(tokCnt==1)
			{
				position_b=atoi(tok);
				//flipping
				//position_e=atoi(tok);
			}
			else if(tokCnt==2)
			{
				position_e=atoi(tok);
				//flipping
				//position_b=atoi(tok);
			}
			else if(tokCnt==3)
			{
				name.append(tok);
			}
			else if(tokCnt==4)
			{
				strand==tok[0];
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		if(chromosome.find("MT")!=std::string::npos)
                {
			continue;
		}
		REGION *regionIn=new REGION;
		regionIn->name=name;
		regionIn->beginning=position_b;
		regionIn->ending=position_e;
		regionIn->chromosome=chromosome;
		regionIn->strand=strand;

		if(regionSet.find(chromosome)==regionSet.end())
		{
			map<string,Framework::REGION*>* in=new map<string,Framework::REGION*>;
			(*in)[name]=regionIn;
			regionSet[chromosome]=in;
		}
		else
		{
			(*regionSet.find(chromosome)->second)[name]=regionIn;
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
	dist=0;
    for(map<string,vector<double>*>::iterator cIter=signal.begin();cIter!=signal.end();cIter++)
    {   
        string chromosome=cIter->first;
        if(chromosomeSizes.find(chromosome)!=chromosomeSizes.end())
        {   
            dist+=chromosomeSizes.find(chromosome)->second;
        }
    }
	double normF=total/dist;
	cout << "normalization by " << normF*scalar << endl;
	for(map<string,vector<double>*>::iterator cIter=signal.begin();cIter!=signal.end();cIter++)
	{
		string chromosome=cIter->first;
		vector<double>* values_fg=signal[cIter->first];
		for(int i=0;i<values_fg->size();i++)
		{
			double normVal=((*values_fg)[i]/normF)*scalar;
			(*values_fg)[i]=normVal;
		}
	}
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
	cout << "Now going over " << regionSet.size() << " chromosomes" << endl;
	for(map<string,map<string,Framework::REGION*>*>::iterator cIter=regionSet.begin();cIter!=regionSet.end();cIter++)
	{
		map<string,Framework::REGION*>* regions=cIter->second;
		cout << "Working on chromosome " << cIter->first << " with " << regions->size() << " regions" << endl;
		if(signal.find(cIter->first)==signal.end())
		{
			cout << "No signal: Skipping chromosome " << cIter->first << endl;
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
	 vector<double>* svals=signal[chromosome];//values for the whole chromosome
		
	//loop over regions
	for(map<string,Framework::REGION*>::iterator rIter=regionSet[chromosome]->begin();rIter!=regionSet[chromosome]->end();rIter++)
	{
		Framework::REGION* r=rIter->second;
		int binstart=r->beginning;
		int binend=r->ending;
		// sanity check against chrom sizes -- if bin goes off the end, make it stop at the last position
		binend=min(binend,chromosomeSizes[chromosome]);
		//print file 
		if (printProfile)
		{
		       sigFile << r->name;
		}

		double mean=0;
		for(int b=binstart;b<=binend;b++)
		{
			int coord_key=b;
			mean=mean+(*svals)[coord_key];
			//print profile
			if (printProfile)
                	{
                       		sigFile << "\t" << (*svals)[coord_key];
                	}	
		}
		if(printProfile)
		{
			sigFile << endl;
		}
		// take avg counts over the motif instance.
		mean=mean/(binend-binstart);
		r->signal=mean;
		avgFile << r->name <<"\t" << r->signal << endl;
	} // end loop over regions for this chromosome;
	return 0;
}

// updated by DC for motifnet, and added profile printing back
// added seq depth norm option
int
main(int argc, const char** argv)
{
	string USAGE="Usage: aggregateSignal bedfile chromosomesizes signal outfilePrefix [-p/--printProfile <windowSize>, -n/--normalize <scalar>]\n"
		"If -p/--printProfile option, will print per-base-pair signal for motif instances. Use with caution as it will create huge files.\n"
		"\tThis requires a parameter for defining a uniformly sized window around each motif center: center +/- window.\n"
		"\tYou need to set the window size at least as wide as half of your largest motif.\n"
		"If -n/--normalize <scalar>, will normalize for sequencing depth and multiple values by any scalar given with this option";

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

	const char* regionFile; //arg 1
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
				case 1: regionFile=argv[a]; break;
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
	int regionsOK=fw.readRegions(regionFile);
	if (regionsOK != 0)
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
	
