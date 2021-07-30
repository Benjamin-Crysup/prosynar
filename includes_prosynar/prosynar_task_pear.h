#ifndef PROSYNAR_TASK_PEAR_H
#define PROSYNAR_TASK_PEAR_H 1

#include "prosynar_task.h"

/**Merge by sliding until the largest number of matches found.*/
class PEARMerger : public ProsynarMerger{
public:
	/**Set up a basic parser.*/
	PEARMerger();
	/**Tear down.*/
	~PEARMerger();
	int handleUnknownArgument(int argc, char** argv, std::ostream* helpOut);
	void printExtraGUIInformation(std::ostream* toPrint);
	int posteriorCheck();
	void initialize(ProsynarArgumentParser* baseArgs);
	int mergePair(int threadInd, CRBSAMFileContents* read1, CRBSAMFileContents* read2, std::string* mergeSeq, std::vector<double>* mergeQual, std::string* errRep);
	
	/**The required number of bases of overlap.*/
	intptr_t reqOverlap;
	/**Whether to reclaim soft clipped bases.*/
	bool softReclaim;
	/**Signifigance level.*/
	double sigLevel;
	/**Whether to merge qualities the original way.*/
	bool origRecipe;
	/**Score gain for a match.*/
	double matchPoints;
	/**Score gain for a mismatch.*/
	double mismatchPoints;
	/**Do the statistical test with the observed overlap.*/
	bool testObsOver;
	
	/**Save the base arguments.*/
	ProsynarArgumentParser* saveArgs;
	/**Places to store the cigar locations.*/
	std::vector< std::vector<intptr_t> > cigLocSet;
	/**Places to store sequences.*/
	std::vector< std::vector<char> > seqASet;
	/**Places to store sequences.*/
	std::vector< std::vector<char> > seqqASet;
	/**Places to store sequences.*/
	std::vector< std::vector<double> > seqqdASet;
	/**Places to store sequences.*/
	std::vector< std::vector<char> > seqBSet;
	/**Places to store sequences.*/
	std::vector< std::vector<char> > seqqBSet;
	/**Places to store sequences.*/
	std::vector< std::vector<double> > seqqdBSet;
	/**Places to store sequences.*/
	std::vector< std::vector<char> > seqBRSet;
	/**Places to store sequences.*/
	std::vector< std::vector<char> > seqqBRSet;
	/**Places to store sequences.*/
	std::vector< std::vector<double> > seqqdBRSet;
	/**Places to store infix.*/
	std::vector< std::vector<int> > seqILSet;
	/**Places to store sequences.*/
	std::vector< std::vector<double> > seqdILSet;
	/**Places to store infix.*/
	std::vector< std::vector<int> > seqIRSet;
	/**Places to store sequences.*/
	std::vector< std::vector<double> > seqdIRSet;
};

/**Factory function.*/
ProsynarMerger* factoryPearMerger();

#endif