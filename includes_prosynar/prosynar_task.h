#ifndef PROSYNAR_TASK_H
#define PROSYNAR_TASK_H 1

#include "whodun_args.h"
#include "whodun_align_affinepd.h"
#include "whodun_parse_table_genome.h"

class ProsynarArgumentParser;

/**Filter out pairs that should not be merged.*/
class ProsynarFilter : public ArgumentParser{
public:
	/**Set up a basic parser.*/
	ProsynarFilter();
	/**Tear down.*/
	virtual ~ProsynarFilter();
	
	/**
	 * Perform initial setup.
	 * @param baseArgs The base arguments and the database.
	 */
	virtual void initialize(ProsynarArgumentParser* baseArgs) = 0;
	/**
	 * Examine two parts of a pair and determine whether a merge should be abandoned.
	 * @param threadInd The thread this is.
	 * @param read1 The first part of the pair.
	 * @param read2 The second part of the pair.
	 * @param errRep The place to put an error message, if any.
	 * @return Whether the two can be merged (1), or should be abandoned (0) or encountered an error (-1).
	 */
	virtual int filterPair(int threadInd, CRBSAMFileContents* read1, CRBSAMFileContents* read2, std::string* errRep) = 0;
};

/**Merger survivors.*/
class ProsynarMerger : public ArgumentParser{
public:
	/**Set up a basic parser.*/
	ProsynarMerger();
	/**Tear down.*/
	virtual ~ProsynarMerger();
	
	/**
	 * Perform initial setup.
	 * @param baseArgs The base arguments and the database.
	 */
	virtual void initialize(ProsynarArgumentParser* baseArgs) = 0;
	/**
	 * Examine two parts of a pair and determine whether a merge should be abandoned.
	 * @param threadInd The thread this is.
	 * @param read1 The first part of the pair.
	 * @param read2 The second part of the pair.
	 * @param mergeSeq The place to put the merged sequence.
	 * @param mergeQual The place to put the merged quality.
	 * @param errRep The place to put an error message, if any.
	 * @return Whether there was an error (-1) or the merge otherwise failed (1).
	 */
	virtual int mergePair(int threadInd, CRBSAMFileContents* read1, CRBSAMFileContents* read2, std::string* mergeSeq, std::vector<double>* mergeQual, std::string* errRep) = 0;
};

/**Parse common arguments for prosynar.*/
class ProsynarArgumentParser : public ArgumentParser{
public:
	/**Set up a basic parser.*/
	ProsynarArgumentParser();
	/**Tear down.*/
	~ProsynarArgumentParser();
	
	/**
	 * Load and let the sub pieces load.
	 */
	void performSetup();
	
	int handleUnknownArgument(int argc, char** argv, std::ostream* helpOut);
	void printExtraGUIInformation(std::ostream* toPrint);
	int posteriorCheck();
	
	/**Standard lock on stderr.*/
	void* errLock;
	/**The file to write merged sequences to.*/
	char* seqOutFile;
	/**The sam file to write merged sequences to.*/
	char* mergeSamOutFile;
	/**The reference file.*/
	char* refFile;
	/**The problem region file.*/
	char* probFile;
	/**The cost specification file.*/
	char* costFile;
	/**The quality mangle file.*/
	char* qualmFile;
	/**The file to write failed sequences to.*/
	char* failDumpFile;
	/**The number of threads to use.*/
	intptr_t numThread;
	/**The names of the sam files to read from.*/
	std::vector<const char*> samNames;
	/**The filters to use.*/
	std::vector<ProsynarFilter*> useFilters;
	/**The merge method to use.*/
	ProsynarMerger* useMerger;
	
	/**the reference sequences: map from reference name to reference sequence*/
	std::map< std::string, std::string > allRefs;
	/**The sam file to write failed merges to.*/
	OutStream* failDumpS;
	/**The sam file to write failed merges to.*/
	TabularWriter* failDumpT;
	/**The sam file to write failed merges to.*/
	CRBSAMFileWriter* failDumpB;
	
	/**The default problematic regions, by reference. Set by filters and mergers*/
	std::map< std::string , std::vector< std::pair<intptr_t,intptr_t> > >* defProbRegMap;
	/**The position dependent cost information.*/
	std::map<std::string,PositionDependentCostKDTree>* defAllRegCosts;
	/**The quality mangle stuff.*/
	std::map<std::string, PositionDependentQualityMangleSet>* defAllQualMangs;
	
	/**The default problematic regions, if provided, by reference.*/
	std::map< std::string , std::vector< std::pair<intptr_t,intptr_t> > > probRegMap;
	/**The default position dependent cost information, if provided.*/
	std::map<std::string,PositionDependentCostKDTree> allRegCosts;
	/**The default quality mangle stuff, if provided.*/
	std::map<std::string, PositionDependentQualityMangleSet> allQualMangs;
	
	/**Storage for the help string.*/
	std::string helpDocStore;
	/**The default output name.*/
	char defOutFN[2];
};

/**
 * Get the names of the prosynar filters, and factories for them.
 * @param toFill The place to put the stuff.
 */
void getAllProsynarFilters(std::map<std::string,ProsynarFilter*(*)()>* toFill);

/**
 * Get the names of the prosynar mergers, and factories for them.
 * @param toFill The place to put the stuff.
 */
void getAllProsynarMergers(std::map<std::string,ProsynarMerger*(*)()>* toFill);

#endif