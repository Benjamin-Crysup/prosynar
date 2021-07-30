#ifndef PROSYNAR_TASK_REFPROB_H
#define PROSYNAR_TASK_REFPROB_H 1

#include "prosynar_task.h"

class ProblematicRegionFilter : public ProsynarFilter{
public:
	/**Set up a default filter.*/
	ProblematicRegionFilter();
	/**Tear down.*/
	~ProblematicRegionFilter();
	int handleUnknownArgument(int argc, char** argv, std::ostream* helpOut);
	void printExtraGUIInformation(std::ostream* toPrint);
	int posteriorCheck();
	void initialize(ProsynarArgumentParser* baseArgs);
	int filterPair(int threadInd, CRBSAMFileContents* read1, CRBSAMFileContents* read2, std::string* errRep);
	
	/**The problem region file.*/
	char* probFile;
	/**The reference cost specification file.*/
	char* costRefFile;
	/**The reference quality mangle file.*/
	char* qualmRefFile;
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
	
	/**Whether the expanded regions are problematic.*/
	std::map<std::string, std::vector<bool> > expandProbRegP;
	/**The regions of intereset.*/
	std::map<std::string, std::vector< std::pair<intptr_t,intptr_t> > > expandProbReg;
	
	/**The problematic regions this uses.*/
	std::map< std::string , std::vector< std::pair<intptr_t,intptr_t> > >* useProbReg;
	/**The position dependent cost information.*/
	std::map<std::string,PositionDependentCostKDTree>* useRegCosts;
	/**The quality mangle stuff.*/
	std::map<std::string, PositionDependentQualityMangleSet>* useQualMangs;
	
	/**The default problematic regions, if provided, by reference.*/
	std::map< std::string , std::vector< std::pair<intptr_t,intptr_t> > > probRegMap;
	/**The position dependent cost information, if provided.*/
	std::map<std::string,PositionDependentCostKDTree> allRegCosts;
	/**The quality mangle stuff, if provided.*/
	std::map<std::string, PositionDependentQualityMangleSet> allQualMangs;
	
	/**Places to store the cigar locations.*/
	std::vector< std::vector<intptr_t> > cigLocSet1;
	/**Places to store the cigar locations.*/
	std::vector< std::vector<intptr_t> > cigLocSet2;
	/**Places to store names.*/
	std::vector<std::string> nameTmpSet;
	/**Places to store read sequence.*/
	std::vector<std::string> seqTmpSet;
	/**Places to store read sequence.*/
	std::vector< std::vector<char> > qualTmpSet;
	/**Places to store reference sequence.*/
	std::vector<std::string> refATmpSet;
	/**Places to store reference sequence.*/
	std::vector<std::string> refBTmpSet;
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
ProsynarFilter* factoryProblematicRegionFilter();


#endif