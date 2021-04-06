#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <map>
#include <string>
#include <vector>
using namespace std;
class Framework
{
	public:
	Framework(bool printProfile, int windowSize);
	~Framework();
	int readTSS(const char*);
	int readSignals(const char*);
	int readInput(const char*);

	// Normalize coutn values by total counts and then multiply by scalar.
	int normalize(double scalar);
	
	// Get signal
	int getSignalInRegion(const char*); //, bool motifnet, bool printPerBP);

	int readChromosomeSizes(const char*);

	private:
	int readData(const char* fName,map<string,vector<double>*>& data,double& );

	// Do processing for one chromosome
	int processChromosome(string chromosome, ofstream& avgFile, ofstream& sigFile); 

	// Helper function to score up all possible edges on a chromosome
	int populateTSSRecords(string chromosome, 
			map<string, map<string, double>* >& storedSignal, 
			map<string, int>& storedNumRegions);
	
	// Helper function to write out results
	int reportTSSRecords(string chromosome, 
		map<string, map<string, double>* >& storedSignal, 
		map<string, int>& storedNumRegions,
		ofstream& avgFile,
		ofstream& sigFile);

	// Yet another helper function. 
	// Prints out per-BP signal around the center of each motif instance
	// for the bestTSS for a motif->gene pair.
	// Pads out to the window size specified in the command line.
	int printInstanceProfiles(ofstream& sigFile, string chromosome, string bestTSS);

	// Cleans up storedSignal data 
	// key: motif->gene pair
	// value: key: TSS, value=score for that TSS.
	int clearStoredSignal(map<string, map<string, double>* >& storedSignal);


	
	// we are actually using this data structure for a triplet of motif/gene/TSS.
	// the best scoring TSS will contribute its score to the edge motif->gene.
	struct MotifTSS
	{
		// This is the specific edge: triplet of motif, gene, TSS.
		string edgeTriplet; // save the motif/gene/TSS combination here
		string motifGenePair; // save the gene-motif combination here.
		map<int,int> tsses; // are we using this?
		char strand;
		
		// paired vectors - each index is a different motif instance.
		vector<int> beginnings;
		vector<int> endings;		
	};


	map<string,map<int,double>*> binInfo;
	
	// we initially flatten out all possible TSSes for each motif->gene edge
	map<string,map<string,MotifTSS*>*> tssSet; //{ chromosome : { Motif->(gene.TSS) : motifTSS triple } }
	map<string,map<int,MotifTSS*>*> tsscoordSet; //{ chromosome : { startPos : motifTSS triple } }

	
	map<string,vector<double>*> signal; // preallocates everything up front

	// legacy
	//map<string,vector<double>*> bg;
	//map<string,vector<double>*> logRatios;
	double dist;
	double total; 	
	double minVal;
	double signalAvg;
	double inputAvg;
	double totalSignal;
	
	map<string,int> chromosomeSizes;

	bool motifNet;
	bool printProfile;
	int windowSize; // defines a window 2x windowSize around the motif center for profile printing
	double normF;
	bool doNormalize;
};
#endif
