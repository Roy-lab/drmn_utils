#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <iostream>
#include <map>
#include <string>
#include <vector>
using namespace std;
class HierarchicalCluster;
class OptimalLeafOrder;
class GeneExpManager;
class MappedOrthogroupReader;
class Framework
{
	public:
	Framework();
	~Framework();
	int readDataMatrix(const char*);
	int readDataMatrixCellTypeSpecific(const char*); // read from the speciesspecnames_clusterassign files, a different dir per celltype
	int readOGIDs(const char*,const char*);
	int setSrcCellType(const char*);
	int readOrder(const char* order, const char* result); // reads celltype order 
	int generateTransitioningGeneSets(double threshold,const char* outdir,int mingenesetSize);
	int generateOrderedClusterMeans(const char* aFName);
	int reorder(const char* aFSuff,double threshold,const char* outdir,int mingenesetSize);
	private:
	int getTree();
	int optimalLeafOrder();
	int readColumns(char*);
	int createMarkHeaders();
	vector<double>* getMeanContinuous(vector<string>*);
	vector<int>* getMeanDiscrete(vector<string>*);
	const char* getNameinCelltype(const char* locus,const char* celltype);

	// map the cell type specific name, locus (in thisCelltype), to the source celltype
	const char* getSourceName(const char* locus, const char* thisCelltype);

	HierarchicalCluster cluster;
	OptimalLeafOrder olo;
	MappedOrthogroupReader mor;
	map<string,HierarchicalClusterNode*> nodeSet;
	map<string,HierarchicalClusterNode*> backup;
	map<string,int> nameIDMap;
	map<string,int> attribNameIDMap;
	map<int,string> attribIDNameMap;
	map<string,int> attribNameIDMap_Profile;
	map<int,string> attribIDNameMap_Profile;
	map<string,GeneExpManager*> markProfileSet;
	map<int,vector<string>*> clusterset;
	string srcCelltype;

	vector<int> hasData; // which celltypes had data
};
#endif
