#ifndef PROSYNAR_TASK_REFOVER_H
#define PROSYNAR_TASK_REFOVER_H 1

#include "prosynar_task.h"

/**Filter on simple reference overlap.*/
class ReferenceOverlapFilter : public ProsynarFilter{
public:
	/**Set up a default filter.*/
	ReferenceOverlapFilter();
	/**Tear down.*/
	~ReferenceOverlapFilter();
	int handleUnknownArgument(int argc, char** argv, std::ostream* helpOut);
	void printExtraGUIInformation(std::ostream* toPrint);
	int posteriorCheck();
	void initialize(ProsynarArgumentParser* baseArgs);
	int filterPair(int threadInd, CRBSAMFileContents* read1, CRBSAMFileContents* read2, std::string* errRep);
	
	/**The required number of bases of overlap.*/
	intptr_t reqOverlap;
	/**Save the base arguments.*/
	ProsynarArgumentParser* saveArgs;
	
	/**Places to store the cigar locations.*/
	std::vector< std::vector<intptr_t> > cigLocSet;
	/**Places to store names.*/
	std::vector<std::string> nameTmpSet;
};

/**Factory function.*/
ProsynarFilter* factoryReferenceOverlapFilter();

class ProbabilisticReferenceOverlapFilter : public ProsynarFilter{
public:
	/**Set up a default filter.*/
	ProbabilisticReferenceOverlapFilter();
	/**Tear down.*/
	~ProbabilisticReferenceOverlapFilter();
	int handleUnknownArgument(int argc, char** argv, std::ostream* helpOut);
	void printExtraGUIInformation(std::ostream* toPrint);
	int posteriorCheck();
	void initialize(ProsynarArgumentParser* baseArgs);
	int filterPair(int threadInd, CRBSAMFileContents* read1, CRBSAMFileContents* read2, std::string* errRep);
	
	/**The reference cost specification file.*/
	char* costRefFile;
	/**The reference quality mangle file.*/
	char* qualmRefFile;
	/**The required number of bases of overlap.*/
	intptr_t reqOverlap;
	/**The threshold likelihood ratio.*/
	double threshLR;
	/**The worst rank to look at.*/
	intptr_t uptoRank;
	/**The number of alignments to report.*/
	intptr_t uptoCount;
	/**The amount of overrun for region realign.*/
	intptr_t overRun;
	/**The number of bases to extend per soft clipped base.*/
	intptr_t softReclaim;
	/**The maximum number of times to see a score before stopping.*/
	intptr_t hotfuzz;
	/**The probability of opening a gap.*/
	double lproGapOpen;
	/**The probability of extending a gap.*/
	double lproGapExtend;
	
	/**Save the base arguments.*/
	ProsynarArgumentParser* saveArgs;
	
	/**The position dependent cost information.*/
	std::map<std::string,PositionDependentCostKDTree>* useRegCosts;
	/**The quality mangle stuff.*/
	std::map<std::string, PositionDependentQualityMangleSet>* useQualMangs;
	
	/**The position dependent cost information, if provided.*/
	std::map<std::string,PositionDependentCostKDTree> allRegCosts;
	/**The quality mangle stuff, if provided.*/
	std::map<std::string, PositionDependentQualityMangleSet> allQualMangs;
	
	/**Places to store names.*/
	std::vector<std::string> nameTmpSet;
	/**Places to store alignment bounds.*/
	std::vector< std::vector< std::pair<intptr_t,intptr_t> > > alnBndsSet1;
	/**Places to store alignment probabilities.*/
	std::vector< std::vector<double> > alnProsSet1;
	/**Places to store alignment bounds.*/
	std::vector< std::vector< std::pair<intptr_t,intptr_t> > > alnBndsSet2;
	/**Places to store alignment probabilities.*/
	std::vector< std::vector<double> > alnProsSet2;
	/**Places to store the cigar locations.*/
	std::vector< std::vector<intptr_t> > cigLocSet;
	/**Places to store reference.*/
	std::vector<std::string> refTmpSet;
	/**Places to store read.*/
	std::vector<std::string> readTmpSet;
	/**Places to store read qvalues.*/
	std::vector< std::vector<char> > readQTmpSet;
	/**Places to store the scores.*/
	std::vector< std::vector<intptr_t> > scoreSet;
	/**Places to store the scores.*/
	std::vector< std::vector<uintptr_t> > scoreSeenSet;
	/**Places to store read qvalues.*/
	std::vector< std::vector<double> > readQPTmpSet;
	/**Save rebased costs.*/
	std::vector<PositionDependentCostKDTree> rebaseCosts;
	/**Save mangled costs.*/
	std::vector<PositionDependentCostKDTree> mangleCosts;
	/**Save rebased mangles.*/
	std::vector<PositionDependentQualityMangleSet> rebaseMangs;
	/**Save alignments.*/
	std::vector<PositionDependentAffineGapLinearPairwiseAlignment> workAlns;
	/**Save iteration tokens.*/
	std::vector<LinearPairwiseAlignmentIteration*> workIters;
};

/**Factory function.*/
ProsynarFilter* factoryProbabilisticReferenceOverlapFilter();

#endif