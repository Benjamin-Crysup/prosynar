#include "prosynar_task_refover.h"

#include <string.h>
#include <algorithm>

#include "whodun_probutil.h"
#include "whodun_parse_seq.h"

ReferenceOverlapFilter::ReferenceOverlapFilter(){
	reqOverlap = 1;
	myMainDoc = "prosynar -- Frover [OPTION]\nFilter by the amount of overlap implied by the initial alignments.\nThe OPTIONS are:\n";
	myVersionDoc = "ProSynAr Frover 1.0";
	myCopyrightDoc = "Copyright (C) 2019 UNT HSC Center for Human Identification";
	ArgumentParserIntMeta overMeta("Overlap Threshold");
		addIntegerOption("--over", &reqOverlap, 0, "    Specify the number of bases of overlap to require.\n    --over 10\n", &overMeta);
}

ReferenceOverlapFilter::~ReferenceOverlapFilter(){
}

int ReferenceOverlapFilter::handleUnknownArgument(int argc, char** argv, std::ostream* helpOut){
	if(strcmp(argv[0],"--")==0){
		return 0;
	}
	argumentError.append("Unknown argument ");
	argumentError.append(argv[0]);
	return -1;
}

void ReferenceOverlapFilter::printExtraGUIInformation(std::ostream* toPrint){
	//nothing special
}

int ReferenceOverlapFilter::posteriorCheck(){
	if(reqOverlap < 0){
		argumentError = "Overlap threshold must be non-negative.";
		return 1;
	}
	return 0;
}

void ReferenceOverlapFilter::initialize(ProsynarArgumentParser* baseArgs){
	saveArgs = baseArgs;
	cigLocSet.resize(baseArgs->numThread);
	nameTmpSet.resize(baseArgs->numThread);
}

int ReferenceOverlapFilter::filterPair(int threadInd, CRBSAMFileContents* read1, CRBSAMFileContents* read2, std::string* errRep){
	std::vector<intptr_t>* cigVec = &(cigLocSet[threadInd]);
	//should both be mapped
	if((read1->entryPos < 0) || (read2->entryPos < 0)){ return 0; }
	if((read1->entryReference.size() == 0) || (read2->entryReference.size() == 0)){ return 0; }
	//should both be on the same reference
	if((read1->entryReference.size() != read2->entryReference.size()) || (memcmp(&(read1->entryReference[0]), &(read2->entryReference[0]), read1->entryReference.size())!=0)){ return 0; }
	//expand the cigars, figure the bounds
	std::pair<intptr_t,intptr_t> read1B;
	try{
		cigVec->clear();
		/*std::pair<uintptr_t,uintptr_t> mainSClip = */cigarStringToReferencePositions(read1->entryPos, &(read1->entryCigar), cigVec);;
		read1B = getCigarReferenceBounds(cigVec);
	}catch(std::exception& err){
		errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
		errRep->push_back(':');
		errRep->push_back(' ');
		errRep->append(err.what());
		return 0;
	}
	std::pair<intptr_t,intptr_t> read2B;
	try{
		cigVec->clear();
		/*std::pair<uintptr_t,uintptr_t> mainSClip = */cigarStringToReferencePositions(read2->entryPos, &(read2->entryCigar), cigVec);;
		read2B = getCigarReferenceBounds(cigVec);
	}catch(std::exception& err){
		errRep->insert(errRep->end(),read2->entryName.begin(), read2->entryName.end());
		errRep->push_back(':');
		errRep->push_back(' ');
		errRep->append(err.what());
		return 0;
	}
	//make sure they are mapped
	if(read1B.first < 0){ return 0; }
	if(read2B.first < 0){ return 0; }
	//check for any overlap
	if(read1B.second < read2B.first){ return 0; }
	if(read1B.first > read2B.second){ return 0; }
	//calculate the total overlap
	intptr_t totOver = std::min(read1B.second, read2B.second) - std::max(read1B.first, read2B.first);
	if(totOver < reqOverlap){
		return 0;
	}
	return 1;
}

ProsynarFilter* factoryReferenceOverlapFilter(){
	return new ReferenceOverlapFilter();
}


ProbabilisticReferenceOverlapFilter::ProbabilisticReferenceOverlapFilter(){
	costRefFile = 0;
	qualmRefFile = 0;
	reqOverlap = 1;
	threshLR = 0.0;
	uptoRank = 4;
	uptoCount = 100;
	lproGapOpen = -3.0;
	lproGapExtend = -1.0;
	overRun = 20;
	softReclaim = 1;
	hotfuzz = 1000;
	myMainDoc = "prosynar -- Fprover [OPTION]\nFilter by the amount of overlap of alignments, weighted by probability.\nThe OPTIONS are:\n";
	myVersionDoc = "ProSynAr Fprover 1.0";
	myCopyrightDoc = "Copyright (C) 2019 UNT HSC Center for Human Identification";
	ArgumentParserIntMeta overMeta("Overlap Threshold");
		addIntegerOption("--over", &reqOverlap, 0, "    Specify the number of bases of overlap to require.\n    --over 10\n", &overMeta);
	ArgumentParserStrMeta costMeta("Alignment Parameter File");
		costMeta.isFile = true;
		costMeta.fileExts.insert(".pdc");
		addStringOption("--cost", &costRefFile, 0, "    Specify the position dependent cost function.\n    --cost File.pdc\n", &costMeta);
	ArgumentParserIntMeta expanMeta("Realign Region Expansion");
		addIntegerOption("--overrun", &overRun, 0, "    Specify the number of extra bases to get when realigning.\n    --overrun 20\n", &expanMeta);
	ArgumentParserIntMeta rankMeta("Non-optimal Scores");
		addIntegerOption("--rank", &uptoRank, 0, "    Specify the number of extra scores to report (0 for best only).\n    --rank 4\n", &rankMeta);
	ArgumentParserIntMeta countMeta("Variant Alignment Count");
		addIntegerOption("--count", &uptoCount, 0, "    Specify the maximum number of alignments to report.\n    Best scores given precedence, zero to report all alternatives.\n    --count 100\n", &countMeta);
	ArgumentParserIntMeta recSoftMeta("Reclaim Soft-clip Compensation");
		addIntegerOption("--reclaimsoft", &softReclaim, 0, "    Use soft clipped bases when realigning, with overrun per base.\n    --reclaimsofft 1\n", &recSoftMeta);
	ArgumentParserStrMeta qualmMeta("Quality Mangle File");
		qualmMeta.isFile = true;
		qualmMeta.fileExts.insert(".qualm");
		addStringOption("--qualm", &qualmRefFile, 0, "    Specify how to modify alignment parameters using their quality.\n    --qualm File.qualm\n", &qualmMeta);
	ArgumentParserIntMeta hfuzzMeta("Score Encounter Limit");
		addIntegerOption("--hfuzz", &hotfuzz, 0, "    The maximum number of times to examine a score during score search.\n    Zero for a FULL traversal.\n    --hfuzz 1000\n", &hfuzzMeta);
	ArgumentParserFltMeta thrLRMeta("Liklihood Threshold");
		addFloatOption("--lr", &threshLR, 0, "    The (log10) liklihood ratio threshold on overlap.\n    --lr 0.0\n", &thrLRMeta);
	ArgumentParserFltMeta gapOMeta("Gap Open Probability");
		addFloatOption("--gopen", &lproGapOpen, 0, "    The (log10) probability of opening a gap.\n    --gopen -3.0\n", &gapOMeta);
	ArgumentParserFltMeta gapEMeta("Gap Extend Probability");
		addFloatOption("--gext", &lproGapExtend, 0, "    The (log10) probability of extending a gap.\n    --gext -1.0\n", &gapEMeta);
}

ProbabilisticReferenceOverlapFilter::~ProbabilisticReferenceOverlapFilter(){
	for(uintptr_t i = 0; i<workIters.size(); i++){
		if(workIters[i]){ delete(workIters[i]); }
	}
	workIters.clear();
}

int ProbabilisticReferenceOverlapFilter::handleUnknownArgument(int argc, char** argv, std::ostream* helpOut){
	if(strcmp(argv[0],"--")==0){
		return 0;
	}
	argumentError.append("Unknown argument ");
	argumentError.append(argv[0]);
	return -1;
}

void ProbabilisticReferenceOverlapFilter::printExtraGUIInformation(std::ostream* toPrint){
	//nothing special
}

int ProbabilisticReferenceOverlapFilter::posteriorCheck(){
	if(reqOverlap < 0){
		argumentError = "Overlap threshold must be non-negative.";
		return 1;
	}
	if(overRun < 0){
		argumentError = "overrun must be non-negative.";
		return 1;
	}
	if(uptoRank < 0){
		argumentError = "rank must be non-negative.";
		return 1;
	}
	if(uptoCount < 0){
		argumentError = "count must be non-negative.";
		return 1;
	}
	if(softReclaim < 0){
		argumentError = "reclaimsofft must be non-negative.";
		return 1;
	}
	if(hotfuzz < 0){
		argumentError = "hfuzz must be non-negative.";
		return 1;
	}
	return 0;
}

void ProbabilisticReferenceOverlapFilter::initialize(ProsynarArgumentParser* baseArgs){
	saveArgs = baseArgs;
	//figure out the costs
	if(costRefFile){
		std::string fileConts;
		FileInStream redFile(costRefFile);
		readStream(&redFile, &fileConts);
		parseMultiregionPositionDependentCost(fileConts.c_str(), fileConts.c_str() + fileConts.size(), &allRegCosts);
		useRegCosts = &allRegCosts;
		if(saveArgs->defAllRegCosts == 0){
			saveArgs->defAllRegCosts = &allRegCosts;
		}
	}
	else if(saveArgs->defAllRegCosts){
		useRegCosts = saveArgs->defAllRegCosts;
	}
	else{
		throw std::runtime_error("No alignment costs specified.");
	}
	//quality mangles
	if(qualmRefFile){
		std::string fileConts;
		FileInStream redFile(qualmRefFile);
		readStream(&redFile, &fileConts);
		parseMultiregionPositionQualityMangle(fileConts.c_str(), fileConts.c_str() + fileConts.size(), &allQualMangs);
		useQualMangs = &allQualMangs;
		if(saveArgs->defAllQualMangs == 0){
			saveArgs->defAllQualMangs = &allQualMangs;
		}
	}
	else if(saveArgs->defAllQualMangs){
		useQualMangs = saveArgs->defAllQualMangs;
	}
	else{
		useQualMangs = &allQualMangs;
	}
	//prepare thread local storage
	nameTmpSet.resize(baseArgs->numThread);
	alnBndsSet1.resize(baseArgs->numThread);
	alnProsSet1.resize(baseArgs->numThread);
	alnBndsSet2.resize(baseArgs->numThread);
	alnProsSet2.resize(baseArgs->numThread);
	cigLocSet.resize(baseArgs->numThread);
	refTmpSet.resize(baseArgs->numThread);
	readTmpSet.resize(baseArgs->numThread);
	readQTmpSet.resize(baseArgs->numThread);
	scoreSet.resize(baseArgs->numThread);
	scoreSeenSet.resize(baseArgs->numThread);
	readQPTmpSet.resize(baseArgs->numThread);
	for(uintptr_t i = 0; i<(uintptr_t)(baseArgs->numThread); i++){
		scoreSet[i].resize(uptoRank+1);
	}
	rebaseCosts.resize(baseArgs->numThread);
	mangleCosts.resize(baseArgs->numThread);
	rebaseMangs.resize(baseArgs->numThread);
	workAlns.resize(baseArgs->numThread);
	LinearPairwiseAlignmentIteration* nullIt = 0;
	workIters.insert(workIters.end(), baseArgs->numThread, nullIt);
}

/**
 * Get the alignments for an entry.
 * @param threadInd The thread index.
 * @param baseFil The base filter: get options.
 * @param read1 The read in question.
 * @param fillBnd The place to put the bounds of the alignments.
 * @param fillPro The place to put the probabilities of the alignments.
 * @param errRep The place to note error messages.
 * @return Whether there was a problem.
 */
int probabilisticROFGetAlignments(int threadInd, ProbabilisticReferenceOverlapFilter* baseFil, CRBSAMFileContents* read1, std::vector< std::pair<intptr_t,intptr_t> >* fillBnd, std::vector<double>* fillPro, std::string* errRep){
	fillBnd->clear();
	fillPro->clear();
	//common values
		uintptr_t maxCount = baseFil->uptoCount;
		intptr_t worstScore = ((intptr_t)-1) << (8*sizeof(intptr_t) - 1);
		std::greater<intptr_t> compMeth;
	//thread local storage
		std::string* nameTmp = &(baseFil->nameTmpSet[threadInd]);
		std::vector<intptr_t>* cigVec = &(baseFil->cigLocSet[threadInd]);
		std::string* mainRefSeq = &(baseFil->refTmpSet[threadInd]);
		std::string* mainReadSeq = &(baseFil->readTmpSet[threadInd]);
		std::vector<char>* mainQualStore = &(baseFil->readQTmpSet[threadInd]);
		std::vector<intptr_t>* mainScores = &(baseFil->scoreSet[threadInd]);
		std::vector<uintptr_t>* packScoreSeen = &(baseFil->scoreSeenSet[threadInd]);
		std::vector<double>* mainQualPStore = &(baseFil->readQPTmpSet[threadInd]);
	//idiot check
		if(read1->entrySeq.size() == 0){
			errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
			errRep->append(" missing sequence.");
			return 1;
		}
		if(read1->entryQual.size() == 0){
			errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
			errRep->append(" missing quality.");
			return 1;
		}
		if(read1->entrySeq.size() != read1->entryQual.size()){
			errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
			errRep->append(" has mismatched sequence and quality.");
			return 1;
		}
	//get the relevant reference
		nameTmp->clear(); nameTmp->insert(nameTmp->end(), read1->entryReference.begin(), read1->entryReference.end());
		std::map<std::string,std::string>::iterator refIt = baseFil->saveArgs->allRefs.find(*nameTmp);
		if(refIt == baseFil->saveArgs->allRefs.end()){
			errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
			errRep->append(" sequence for reference ");
			errRep->insert(errRep->end(),read1->entryReference.begin(), read1->entryReference.end());
			errRep->append("not found.");
			return 1;
		}
		std::string* selRef = &(refIt->second);
	//get the cost info
		std::map<std::string,PositionDependentCostKDTree>::iterator costIt = baseFil->useRegCosts->find(*nameTmp);
		if(costIt == baseFil->useRegCosts->end()){
			errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
			errRep->append(" alignment parameters for reference ");
			errRep->insert(errRep->end(),read1->entryReference.begin(), read1->entryReference.end());
			errRep->append("not found.");
			return 1;
		}
		PositionDependentCostKDTree* selCost = &(costIt->second);
	//get the quality mangle
		std::map<std::string,PositionDependentQualityMangleSet>::iterator qualmIt = baseFil->useQualMangs->find(*nameTmp);
		PositionDependentQualityMangleSet* selMang = 0;
		if(qualmIt != baseFil->useQualMangs->end()){
			selMang = &(qualmIt->second);
		}
	//get the bounds of the thing
		std::pair<intptr_t,intptr_t> read1B;
		std::pair<uintptr_t,uintptr_t> mainSClip;
		try{
			cigVec->clear();
			mainSClip = cigarStringToReferencePositions(read1->entryPos, &(read1->entryCigar), cigVec);
			read1B = getCigarReferenceBounds(cigVec);
		}catch(std::exception& err){
			errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
			errRep->push_back(':');
			errRep->push_back(' ');
			errRep->append(err.what());
			return 0;
		}
	//require a map
		if(read1B.first < 0){ return 0; }
	//expand
		read1B.first = std::max((intptr_t)0, read1B.first - (intptr_t)(baseFil->overRun + baseFil->softReclaim*mainSClip.first));
		read1B.second = std::min(selRef->size()-1, read1B.second + (baseFil->overRun + baseFil->softReclaim*mainSClip.second));
		read1B.second++;
	//get the quality
		mainQualStore->clear(); mainQualStore->insert(mainQualStore->end(), read1->entryQual.begin() + (baseFil->softReclaim ? 0 : mainSClip.first), read1->entryQual.end() - (baseFil->softReclaim ? 0 : mainSClip.second));
		mainQualPStore->resize(mainQualStore->size());
		fastaPhredsToLog10Prob(mainQualStore->size(), (unsigned char*)(&((*mainQualStore)[0])), &((*mainQualPStore)[0]));
	//prepare an alignment
		mainRefSeq->clear(); mainRefSeq->insert(mainRefSeq->end(), selRef->begin() + read1B.first, selRef->begin() + read1B.second);
		mainReadSeq->clear(); mainReadSeq->insert(mainReadSeq->end(), read1->entrySeq.begin() + (baseFil->softReclaim ? 0 : mainSClip.first), read1->entrySeq.end() - (baseFil->softReclaim ? 0 : mainSClip.second));
		PositionDependentCostKDTree* mainUseCost = &(baseFil->rebaseCosts[threadInd]);
			mainUseCost->regionsRebased(selCost, read1B.first, read1B.second, -1, -1);
		if(selMang){
			PositionDependentQualityMangleSet* subMang = &(baseFil->rebaseMangs[threadInd]);
				subMang->rebase(selMang, read1B.first, read1B.second);
			PositionDependentCostKDTree* tmpUse = &(baseFil->mangleCosts[threadInd]);
				tmpUse->regionsQualityMangled(mainUseCost, subMang, mainQualStore);
			mainUseCost = tmpUse;
		}
		mainUseCost->produceFromRegions();
		PositionDependentAffineGapLinearPairwiseAlignment* mainAln = &(baseFil->workAlns[threadInd]);
		mainAln->changeProblem(2, mainRefSeq, mainReadSeq, mainUseCost);
		mainAln->prepareAlignmentStructure();
		LinearPairwiseAlignmentIteration* mainIter = baseFil->workIters[threadInd];
		if(!mainIter){
			baseFil->workIters[threadInd] = mainAln->getIteratorToken();
			mainIter = baseFil->workIters[threadInd];
		}
	//get the relevant scores
		int mainNumScore = mainAln->findAlignmentScores(mainIter, mainScores->size(), &((*mainScores)[0]), worstScore, baseFil->hotfuzz);
	//iterate
		packScoreSeen->resize(mainNumScore);
		for(int i = 0; i<mainNumScore; i++){ (*packScoreSeen)[i] = 0; }
		mainAln->startFuzzyIteration(mainIter, (*mainScores)[mainNumScore-1], baseFil->hotfuzz, mainNumScore);
		while(true){
			if(!(mainIter->getNextAlignment())){ break; }
			uintptr_t scoreInd = std::lower_bound(mainScores->begin(), mainScores->begin() + mainNumScore, mainIter->alnScore, compMeth) - mainScores->begin();
			if(maxCount && ((*packScoreSeen)[scoreInd] >= maxCount)){ continue; }
			double curLPro = linearReferenceAlignProbabilityAffine(mainAln, &((*mainQualPStore)[0]), mainIter, baseFil->lproGapOpen, baseFil->lproGapExtend);
			intptr_t lowInd = mainIter->aInds[0] + read1B.first;
			intptr_t higInd = mainIter->aInds[mainIter->aInds.size()-1] + read1B.first;
			fillBnd->push_back(std::pair<intptr_t,intptr_t>(lowInd,higInd));
			fillPro->push_back(curLPro);
			(*packScoreSeen)[scoreInd]++;
			if((scoreInd == 0) && maxCount && ((*packScoreSeen)[scoreInd] >= maxCount)){ break; }
		}
	return 0;
}

int ProbabilisticReferenceOverlapFilter::filterPair(int threadInd, CRBSAMFileContents* read1, CRBSAMFileContents* read2, std::string* errRep){
	//should both be mapped
	if((read1->entryPos < 0) || (read2->entryPos < 0)){ return 0; }
	if((read1->entryReference.size() == 0) || (read2->entryReference.size() == 0)){ return 0; }
	//should both be on the same reference
	if((read1->entryReference.size() != read2->entryReference.size()) || (memcmp(&(read1->entryReference[0]), &(read2->entryReference[0]), read1->entryReference.size())!=0)){ return 0; }
	//get the alignments for both
	std::vector< std::pair<intptr_t,intptr_t> >* alnBnd1 = &(alnBndsSet1[threadInd]);
	std::vector< std::pair<intptr_t,intptr_t> >* alnBnd2 = &(alnBndsSet2[threadInd]);
	std::vector<double>* alnPro1 = &(alnProsSet1[threadInd]);
	std::vector<double>* alnPro2 = &(alnProsSet2[threadInd]);
	alnBnd1->clear(); alnPro1->clear();
	alnBnd2->clear(); alnPro2->clear();
	if(probabilisticROFGetAlignments(threadInd, this, read1, alnBnd1, alnPro1, errRep)){ return -1; }
	if(probabilisticROFGetAlignments(threadInd, this, read2, alnBnd2, alnPro2, errRep)){ return -1; }
	//calculate the probabilities of the two possibilities
	intptr_t numOver = 0;
	ProbabilitySummation overLike;
	intptr_t numDisj = 0;
	ProbabilitySummation disjLike;
	for(unsigned i = 0; i<alnBnd1->size(); i++){
		for(unsigned j = 0; j<alnBnd2->size(); j++){
			double curProb = (*alnPro1)[i] + (*alnPro2)[j];
			intptr_t overMin = std::max((*alnBnd1)[i].first, (*alnBnd2)[j].first);
			intptr_t overMax = std::min((*alnBnd1)[i].second, (*alnBnd2)[j].second);
			intptr_t over = overMax - overMin;
			if(over >= reqOverlap){
				overLike.addNextLogProb(curProb);
				numOver++;
			}
			else{
				disjLike.addNextLogProb(curProb);
				numDisj++;
			}
		}
	}
	//decide what to do
	if(numOver == 0){ return 0; }
	if(numDisj == 0){ return 1; }
	if((overLike.getFinalLogSum() - disjLike.getFinalLogSum()) >= threshLR){
		return 1;
	}
	return 0;
}

ProsynarFilter* factoryProbabilisticReferenceOverlapFilter(){
	return new ProbabilisticReferenceOverlapFilter();
}

