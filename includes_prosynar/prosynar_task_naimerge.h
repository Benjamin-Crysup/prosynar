#ifndef PROSYNAR_TASK_NAIMERGE_H
#define PROSYNAR_TASK_NAIMERGE_H 1

#include "prosynar_task.h"

/**
 * Merge two sequences with expanded indels (-1 in the sequence): both sequences have basewise correspondence.
 * @param seqA THe first sequence.
 * @param qualA The qualities.
 * @param seqB THe second sequence.
 * @param qualB The qualities.
 * @param mergeSeq The place to put the merged sequence.
 * @param mergeQual The place to put the merged quality.
 */
void mergeIndelExpandedSequences(std::vector<int>* seqA, std::vector<double>* qualA, std::vector<int>* seqB, std::vector<double>* qualB, std::string* mergeSeq, std::vector<double>* mergeQual);

/**Merge by alignment, with a simple overlap threshold.*/
class SimpleAlignMerger : public ProsynarMerger{
public:
	/**Set up a basic parser.*/
	SimpleAlignMerger();
	/**Tear down.*/
	~SimpleAlignMerger();
	int handleUnknownArgument(int argc, char** argv, std::ostream* helpOut);
	void printExtraGUIInformation(std::ostream* toPrint);
	int posteriorCheck();
	void initialize(ProsynarArgumentParser* baseArgs);
	int mergePair(int threadInd, CRBSAMFileContents* read1, CRBSAMFileContents* read2, std::string* mergeSeq, std::vector<double>* mergeQual, std::string* errRep);
	
	/**The required number of bases of overlap.*/
	intptr_t reqOverlap;
	/**The affine gap cost file.*/
	char* costReadFile;
	/**Whether to reclaim soft clipped bases.*/
	bool softReclaim;
	/**The biquality mangle.*/
	char* qualmFile;
	
	/**Save the base arguments.*/
	ProsynarArgumentParser* saveArgs;
	/**The read cost information.*/
	AlignCostAffine mergeCost;
	/**The merge cost.*/
	PositionDependentCostKDTree bigCost;
	/**The quality mangle.*/
	PositionalBiQualityMangleSet bigMang;
	
	/**Places to store the cigar locations.*/
	std::vector< std::vector<intptr_t> > cigLocSet;
	/**Places to store mangled costs.*/
	std::vector<PositionDependentCostKDTree> costSet;
	/**Places to store sequences.*/
	std::vector< std::string > seqASet;
	/**Places to store sequences.*/
	std::vector< std::vector<char> > seqqASet;
	/**Places to store sequences.*/
	std::vector< std::vector<double> > seqqdASet;
	/**Places to store sequences.*/
	std::vector< std::string > seqBSet;
	/**Places to store sequences.*/
	std::vector< std::vector<char> > seqqBSet;
	/**Places to store sequences.*/
	std::vector< std::vector<double> > seqqdBSet;
	/**Places to store infix.*/
	std::vector< std::vector<int> > seqILSet;
	/**Places to store sequences.*/
	std::vector< std::vector<double> > seqdILSet;
	/**Places to store infix.*/
	std::vector< std::vector<int> > seqIRSet;
	/**Places to store sequences.*/
	std::vector< std::vector<double> > seqdIRSet;
	/**Saved iterations.*/
	std::vector<LinearPairwiseAlignmentIteration*> runIters;
	/**Place to store running alignments.*/
	std::vector<PositionDependentAffineGapLinearPairwiseAlignment> saveAlns;
};

/**Factory function.*/
ProsynarMerger* factorySimpleAlignMerger();

#endif