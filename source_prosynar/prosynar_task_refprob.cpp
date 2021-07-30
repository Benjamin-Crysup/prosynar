#include "prosynar_task_refprob.h"

#include <string.h>
#include <algorithm>

#include "whodun_probutil.h"
#include "whodun_parse_seq.h"

ProblematicRegionFilter::ProblematicRegionFilter(){
	probFile = 0;
	costRefFile = 0;
	qualmRefFile = 0;
	threshLR = 0.0;
	uptoRank = 4;
	uptoCount = 100;
	lproGapOpen = -3.0;
	lproGapExtend = -1.0;
	overRun = 20;
	softReclaim = 1;
	hotfuzz = 1000;
	myMainDoc = "prosynar -- Fpprobreg [OPTION]\nFilter by the amount of overlap of alignments, weighted by probability.\nThe OPTIONS are:\n";
	myVersionDoc = "ProSynAr Fpprobreg 1.0";
	myCopyrightDoc = "Copyright (C) 2019 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta probMeta("Problematic Region File");
		probMeta.isFile = true;
		probMeta.fileExts.insert(".bed");
		addStringOption("--prob", &probFile, 0, "    Specify the problematic regions in a bed file.\n    --prob File.bed\n", &probMeta);
	ArgumentParserStrMeta costMeta("Alignment Parameter File");
		costMeta.isFile = true;
		costMeta.fileExts.insert(".pdc");
		addStringOption("--cost", &costRefFile, 0, "    Specify the position dependent cost function.\n    --cost File.pdc\n", &costMeta);
	ArgumentParserIntMeta overMeta("Realign Region Expansion");
		addIntegerOption("--overrun", &overRun, 0, "    Specify the number of extra bases to get when realigning.\n    --overrun 20\n", &overMeta);
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

ProblematicRegionFilter::~ProblematicRegionFilter(){
	for(uintptr_t i = 0; i<workIters.size(); i++){
		if(workIters[i]){ delete(workIters[i]); }
	}
	workIters.clear();
}

int ProblematicRegionFilter::handleUnknownArgument(int argc, char** argv, std::ostream* helpOut){
	if(strcmp(argv[0],"--")==0){
		return 0;
	}
	argumentError.append("Unknown argument ");
	argumentError.append(argv[0]);
	return -1;
}

void ProblematicRegionFilter::printExtraGUIInformation(std::ostream* toPrint){
	//nothing special
}

int ProblematicRegionFilter::posteriorCheck(){
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

void ProblematicRegionFilter::initialize(ProsynarArgumentParser* baseArgs){
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
	//problem regions
	if(probFile){
		FileInStream redFile(probFile);
		TSVTabularReader tabRead(0, &redFile);
		BedFileReader proBed(&tabRead);
		//make easy to search
		std::map< std::string , std::vector< std::pair<intptr_t,intptr_t> > >::iterator probMIt;
		std::vector< std::pair<intptr_t,intptr_t> >::iterator proRanItA;
		std::vector< std::pair<intptr_t,intptr_t> >::iterator proRanItB;
		for(uintptr_t i = 0; i<proBed.chromosomes.size(); i++){
			std::pair<intptr_t,intptr_t> curRange(proBed.locStarts[i], proBed.locEnds[i]);
			probRegMap[proBed.chromosomes[i]].push_back(curRange);
		}
		//sort the vectors
		for(probMIt = probRegMap.begin(); probMIt != probRegMap.end(); probMIt++){
			std::sort(probMIt->second.begin(), probMIt->second.end());
		}
		//save for default
		useProbReg = &probRegMap;
		if(saveArgs->defProbRegMap == 0){
			saveArgs->defProbRegMap = &probRegMap;
		}
	}
	else if(saveArgs->defProbRegMap){
		useProbReg = saveArgs->defProbRegMap;
	}
	else{
		throw std::runtime_error("No problem regions specified.");
	}
	//expand the problem region info
	for(std::map< std::string , std::vector< std::pair<intptr_t,intptr_t> > >::iterator probIt = useProbReg->begin(); probIt != useProbReg->end(); probIt++){
		std::vector< std::pair<intptr_t,intptr_t> >* curReg = &(probIt->second);
		std::vector<bool>* curExpRegP = &(expandProbRegP[probIt->first]);
		std::vector< std::pair<intptr_t,intptr_t> >* curExpReg = &(expandProbReg[probIt->first]);
		intptr_t prevEnd = 0;
		for(uintptr_t i = 0; i<curReg->size(); i++){
			if(prevEnd != (*curReg)[i].first){
				curExpRegP->push_back(false);
				curExpReg->push_back( std::pair<intptr_t,intptr_t>(prevEnd, (*curReg)[i].first) );
			}
			curExpRegP->push_back(true);
			curExpReg->push_back( (*curReg)[i] );
			prevEnd = (*curReg)[i].second;
		}
		//had better have the reference
		std::map<std::string,std::string>::iterator refIt = saveArgs->allRefs.find(probIt->first);
		if(refIt == saveArgs->allRefs.end()){
			throw new std::runtime_error("Reference sequence " + probIt->first + " no found");
		}
		if((uintptr_t)prevEnd > refIt->second.size()){
			throw new std::runtime_error("Problem region beyond the end of the reference sequence " + probIt->first);
		}
		if((uintptr_t)prevEnd != refIt->second.size()){
			curExpRegP->push_back(false);
			curExpReg->push_back( std::pair<intptr_t,intptr_t>(prevEnd, refIt->second.size()) );
		}
	}
	//prepare thread local storage
	cigLocSet1.resize(baseArgs->numThread);
	cigLocSet2.resize(baseArgs->numThread);
	nameTmpSet.resize(baseArgs->numThread);
	seqTmpSet.resize(baseArgs->numThread);
	qualTmpSet.resize(baseArgs->numThread);
	refATmpSet.resize(baseArgs->numThread);
	refBTmpSet.resize(baseArgs->numThread);
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
 * Compare two pairs by the first element only.
 * @param itemA A pair to compare.
 * @param itemB A pair to compare.
 * @param Whether itemA is less than itemB.
 */
bool pairCompOnlyFirst(const std::pair<intptr_t,intptr_t>& itemA, const std::pair<intptr_t,intptr_t>& itemB){
	return itemA.first < itemB.first;
}

/**
 * Compare two pairs by the second element only.
 * @param itemA A pair to compare.
 * @param itemB A pair to compare.
 * @param Whether itemA is less than itemB.
 */
bool pairCompOnlySecond(const std::pair<intptr_t,intptr_t>& itemA, const std::pair<intptr_t,intptr_t>& itemB){
	return itemA.second < itemB.second;
}

/**
 * Get the likelihood that one sequence came from another.
 * @param alnPro The sequences.
 * @param readQuals The read quality along the read.
 * @param forAln The alignment iteration.
 * @param lproGapOpen The (log10) probability of opening a gap.
 * @param lproGapExtend The (log10) probability of extending a gap.
 * @param numScores The number of scores in play.
 * @param hotScores The scores that will be seen in the alignment.
 * @param scoreConsider The number of times to consider each score. If the first hits zero, will stop. If null, will never stop or skip.
 * @return log_10(p(read|reference))
 */
double linearReferenceSourceProbabilityAffine(LinearPairwiseSequenceAlignment* alnPro, double* readQuals, LinearPairwiseAlignmentIteration* forAln, double lproGapOpen, double lproGapExtend, int numScores, intptr_t* hotScores, uintptr_t* scoreConsider){
	std::greater<intptr_t> compMeth;
	ProbabilitySummation endLike;
	while(forAln->getNextAlignment()){
		uintptr_t scoreInd = 0;
		if(scoreConsider){
			scoreInd = std::lower_bound(hotScores, hotScores + numScores, forAln->alnScore, compMeth) - hotScores;
			if(scoreConsider[scoreInd] == 0){ continue; }
		}
		double curLPro = linearReferenceAlignProbabilityAffine(alnPro, readQuals, forAln, lproGapOpen, lproGapExtend);
		endLike.addNextLogProb(curLPro);
		if(scoreConsider){
			scoreConsider[scoreInd]--;
			if((scoreInd == 0) && (scoreConsider[scoreInd] == 0)){ break; }
		}
	}
	return endLike.getFinalLogSum();
}

int ProblematicRegionFilter::filterPair(int threadInd, CRBSAMFileContents* read1, CRBSAMFileContents* read2, std::string* errRep){
	intptr_t worstScore = ((intptr_t)-1) << (8*sizeof(intptr_t) - 1);
	std::string* nameTmp = &(nameTmpSet[threadInd]);
	std::vector<intptr_t>* cigVec1 = &(cigLocSet1[threadInd]);
	std::vector<intptr_t>* cigVec2 = &(cigLocSet2[threadInd]);
	std::string* seqTmp = &(seqTmpSet[threadInd]);
	std::vector<char>* qualTmp = &(qualTmpSet[threadInd]);
	std::string* refATmp = &(refATmpSet[threadInd]);
	std::string* refBTmp = &(refBTmpSet[threadInd]);
	std::vector<intptr_t>* mainScores = &(scoreSet[threadInd]);
	std::vector<uintptr_t>* packScoreSeen = &(scoreSeenSet[threadInd]);
	std::vector<double>* mainQualPStore = &(readQPTmpSet[threadInd]);
	//should both be mapped
		if((read1->entryPos < 0) || (read2->entryPos < 0)){ return 0; }
		if((read1->entryReference.size() == 0) || (read2->entryReference.size() == 0)){ return 0; }
	//should both be on the same reference
		if((read1->entryReference.size() != read2->entryReference.size()) || (memcmp(&(read1->entryReference[0]), &(read2->entryReference[0]), read1->entryReference.size())!=0)){ return 0; }
	//idiot check
		if(read1->entrySeq.size() == 0){
			errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
			errRep->append(" missing sequence.");
			return -1;
		}
		if(read1->entryQual.size() == 0){
			errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
			errRep->append(" missing quality.");
			return -1;
		}
		if(read1->entrySeq.size() != read1->entryQual.size()){
			errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
			errRep->append(" has mismatched sequence and quality.");
			return -1;
		}
		if(read2->entrySeq.size() == 0){
			errRep->insert(errRep->end(),read2->entryName.begin(), read2->entryName.end());
			errRep->append(" missing sequence.");
			return -1;
		}
		if(read2->entryQual.size() == 0){
			errRep->insert(errRep->end(),read2->entryName.begin(), read2->entryName.end());
			errRep->append(" missing quality.");
			return -1;
		}
		if(read2->entrySeq.size() != read2->entryQual.size()){
			errRep->insert(errRep->end(),read2->entryName.begin(), read2->entryName.end());
			errRep->append(" has mismatched sequence and quality.");
			return -1;
		}
	//expand the cigars, figure the bounds
		std::pair<intptr_t,intptr_t> read1B;
		std::pair<uintptr_t,uintptr_t> read1SClip;
		try{
			cigVec1->clear();
			read1SClip = cigarStringToReferencePositions(read1->entryPos, &(read1->entryCigar), cigVec1);;
			read1B = getCigarReferenceBounds(cigVec1);
		}catch(std::exception& err){
			errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
			errRep->push_back(':');
			errRep->push_back(' ');
			errRep->append(err.what());
			return -1;
		}
		std::pair<intptr_t,intptr_t> read2B;
		std::pair<uintptr_t,uintptr_t> read2SClip;
		try{
			cigVec2->clear();
			read2SClip = cigarStringToReferencePositions(read2->entryPos, &(read2->entryCigar), cigVec2);;
			read2B = getCigarReferenceBounds(cigVec2);
		}catch(std::exception& err){
			errRep->insert(errRep->end(),read2->entryName.begin(), read2->entryName.end());
			errRep->push_back(':');
			errRep->push_back(' ');
			errRep->append(err.what());
			return -1;
		}
		read1B.second++;
		read2B.second++;
	//make sure they actually are mapped
		if(read1B.first < 0){ return 0; }
		if(read2B.first < 0){ return 0; }
	//swap if out of order
		if((read1B.first > read2B.first) || ((read1B.first == read2B.first) && (read1B.second > read2B.second))){
			return filterPair(threadInd, read2, read1, errRep);
		}
	//get the reference, costs, mangles, and region set
		nameTmp->clear(); nameTmp->insert(nameTmp->end(), read1->entryReference.begin(), read1->entryReference.end());
		std::map<std::string,std::string>::iterator refIt = saveArgs->allRefs.find(*nameTmp);
			if(refIt == saveArgs->allRefs.end()){
				errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
				errRep->append(" sequence for reference ");
				errRep->insert(errRep->end(),read1->entryReference.begin(), read1->entryReference.end());
				errRep->append(" not found.");
				return -1;
			}
			std::string* selRef = &(refIt->second);
		std::map< std::string, std::vector<bool> >::iterator probIt = expandProbRegP.find(*nameTmp);
			if(probIt == expandProbRegP.end()){
				//no problem = no problem
				return 1;
			}
			std::vector<bool>* selProbP = &(probIt->second);
			std::vector< std::pair<intptr_t,intptr_t> >* selProbR = &(expandProbReg[*nameTmp]);
		std::map<std::string,PositionDependentCostKDTree>::iterator costIt = useRegCosts->find(*nameTmp);
			if(costIt == useRegCosts->end()){
				errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
				errRep->append(" alignment parameters for reference ");
				errRep->insert(errRep->end(),read1->entryReference.begin(), read1->entryReference.end());
				errRep->append(" not found.");
				return -1;
			}
			PositionDependentCostKDTree* selCost = &(costIt->second);
		std::map<std::string,PositionDependentQualityMangleSet>::iterator qualmIt = useQualMangs->find(*nameTmp);
			PositionDependentQualityMangleSet* selMang = 0;
			if(qualmIt != useQualMangs->end()){
				selMang = &(qualmIt->second);
			}
	//idiot check the maps to the reference
		if((uintptr_t)(read1B.second) > selRef->size()){
			errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
			errRep->append(" alignment extends beyond reference ");
			errRep->insert(errRep->end(),read1->entryReference.begin(), read1->entryReference.end());
			return -1;
		}
		if((uintptr_t)(read2B.second) > selRef->size()){
			errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
			errRep->append(" alignment extends beyond reference ");
			errRep->insert(errRep->end(),read1->entryReference.begin(), read1->entryReference.end());
			return -1;
		}
	//note which problem regions they are in
		uintptr_t read1LPI = std::upper_bound(selProbR->begin(), selProbR->end(), read1B, pairCompOnlyFirst) - selProbR->begin(); read1LPI--;
		uintptr_t read1HPI = std::lower_bound(selProbR->begin(), selProbR->end(), read1B, pairCompOnlySecond) - selProbR->begin();
		uintptr_t read2LPI = std::upper_bound(selProbR->begin(), selProbR->end(), read2B, pairCompOnlyFirst) - selProbR->begin(); read2LPI--;
		uintptr_t read2HPI = std::lower_bound(selProbR->begin(), selProbR->end(), read2B, pairCompOnlySecond) - selProbR->begin();
	//simple checks first
		switch(read1HPI - read2LPI){
			case 0:
				if((*selProbP)[read1HPI]){
					//start in same bad spot; known bad
					return 0;
				}
				//have to check for cross over
				break;
			case 1:
				//have to check for cross over
				break;
			case 2:
				if(!((*selProbP)[read1HPI-1])){
					//moving both puts them non-repetitive: known good
					return 1;
				}
				//have to check
				break;
			default:
				//big separation, known good
				return 1;
		}
	//check moving the left edge of the right read
		bool canMoveLeftOfRight = (read2LPI != read2HPI);
		if(canMoveLeftOfRight){
			//get the read sequence
				intptr_t breakInd = (*selProbR)[read2LPI].second;
				seqTmp->clear(); qualTmp->clear();
				if(softReclaim){
					seqTmp->insert(seqTmp->end(), read2->entrySeq.begin(), read2->entrySeq.begin() + read2SClip.first);
					qualTmp->insert(qualTmp->end(), read2->entryQual.begin(), read2->entryQual.begin() + read2SClip.first);
				}
				uintptr_t readAlnSize = overRun + softReclaim * read2SClip.first;
				for(uintptr_t i = 0; i<cigVec2->size(); i++){
					if((*cigVec2)[i] >= breakInd){ break; }
					seqTmp->push_back(read2->entrySeq[i + read2SClip.first]);
					qualTmp->push_back(read2->entryQual[i + read2SClip.first]);
					readAlnSize++;
				}
				mainQualPStore->resize(qualTmp->size());
				fastaPhredsToLog10Prob(qualTmp->size(), (unsigned char*)(&((*qualTmp)[0])), &((*mainQualPStore)[0]));
			//get the competing references
				std::pair<uintptr_t,uintptr_t> refAGot(std::max((intptr_t)0, breakInd-(intptr_t)readAlnSize), breakInd);
				std::pair<uintptr_t,uintptr_t> refBGot(breakInd, std::min(selRef->size(), breakInd+readAlnSize));
				refATmp->clear(); refATmp->insert(refATmp->end(), selRef->begin() + refAGot.first, selRef->begin() + refAGot.second);
				refBTmp->clear(); refBTmp->insert(refBTmp->end(), selRef->begin() + refBGot.first, selRef->begin() + refBGot.second);
			//likelihood of A
				PositionDependentCostKDTree* refACosts = &(rebaseCosts[threadInd]);
					refACosts->regionsRebased(selCost, refAGot.first, refAGot.second, -1, -1);
				if(selMang){
					PositionDependentQualityMangleSet* subMang = &(rebaseMangs[threadInd]);
						subMang->rebase(selMang, refAGot.first, refAGot.second);
					PositionDependentCostKDTree* tmpUse = &(mangleCosts[threadInd]);
						tmpUse->regionsQualityMangled(refACosts, subMang, qualTmp);
					refACosts = tmpUse;
				}
				refACosts->produceFromRegions();
				PositionDependentAffineGapLinearPairwiseAlignment* refAAln = &(workAlns[threadInd]);
				refAAln->changeProblem(2, refATmp, seqTmp, refACosts);
				refAAln->prepareAlignmentStructure();
				LinearPairwiseAlignmentIteration* refAIter = workIters[threadInd];
				if(!refAIter){
					workIters[threadInd] = refAAln->getIteratorToken();
					refAIter = workIters[threadInd];
				}
				int numScoreA = refAAln->findAlignmentScores(refAIter, mainScores->size(), &((*mainScores)[0]), worstScore, hotfuzz);
				packScoreSeen->resize(numScoreA);
				for(int i = 0; i<numScoreA; i++){ (*packScoreSeen)[i] = uptoCount; }
				refAAln->startFuzzyIteration(refAIter, (*mainScores)[numScoreA-1], hotfuzz, numScoreA);
				double refALPro = linearReferenceSourceProbabilityAffine(refAAln, &((*mainQualPStore)[0]), refAIter, lproGapOpen, lproGapExtend, numScoreA, &((*mainScores)[0]), uptoCount ? &((*packScoreSeen)[0]) : 0);
			//likelihood of B
				PositionDependentCostKDTree* refBCosts = &(rebaseCosts[threadInd]);
					refBCosts->regionsRebased(selCost, refBGot.first, refBGot.second, -1, -1);
				if(selMang){
					PositionDependentQualityMangleSet* subMang = &(rebaseMangs[threadInd]);
						subMang->rebase(selMang, refBGot.first, refBGot.second);
					PositionDependentCostKDTree* tmpUse = &(mangleCosts[threadInd]);
						tmpUse->regionsQualityMangled(refBCosts, subMang, qualTmp);
					refBCosts = tmpUse;
				}
				refBCosts->produceFromRegions();
				PositionDependentAffineGapLinearPairwiseAlignment* refBAln = &(workAlns[threadInd]);
				refBAln->changeProblem(2, refBTmp, seqTmp, refBCosts);
				refBAln->prepareAlignmentStructure();
				LinearPairwiseAlignmentIteration* refBIter = workIters[threadInd];
				if(!refBIter){
					workIters[threadInd] = refBAln->getIteratorToken();
					refBIter = workIters[threadInd];
				}
				int numScoreB = refBAln->findAlignmentScores(refBIter, mainScores->size(), &((*mainScores)[0]), worstScore, hotfuzz);
				packScoreSeen->resize(numScoreB);
				for(int i = 0; i<numScoreB; i++){ (*packScoreSeen)[i] = uptoCount; }
				refBAln->startFuzzyIteration(refBIter, (*mainScores)[numScoreB-1], hotfuzz, numScoreB);
				double refBLPro = linearReferenceSourceProbabilityAffine(refBAln, &((*mainQualPStore)[0]), refBIter, lproGapOpen, lproGapExtend, numScoreB, &((*mainScores)[0]), uptoCount ? &((*packScoreSeen)[0]) : 0);
			//make the decision
				canMoveLeftOfRight = (refALPro - refBLPro) < threshLR;
		}
	//check moving the right edge of the left read
		bool canMoveRightOfLeft = (read1LPI != read1HPI);
		if(canMoveRightOfLeft){
			//get the read sequence
				intptr_t breakInd = (*selProbR)[read1HPI].first;
				seqTmp->clear(); qualTmp->clear();
				uintptr_t readAlnSize = overRun + softReclaim * read1SClip.second;
				uintptr_t i = cigVec1->size();
				while(i){
					i--;
					intptr_t ccval = (*cigVec1)[i];
					if((ccval >= 0) && (ccval < breakInd)){ break; }
					seqTmp->push_back(read1->entrySeq[i + read1SClip.first]);
					qualTmp->push_back(read1->entryQual[i + read1SClip.first]);
					readAlnSize++;
				}
				std::reverse(seqTmp->begin(), seqTmp->end());
				std::reverse(qualTmp->begin(), qualTmp->end());
				if(softReclaim){
					seqTmp->insert(seqTmp->end(), read1->entrySeq.end()-read1SClip.second, read1->entrySeq.end());
					qualTmp->insert(qualTmp->end(), read1->entryQual.end()-read1SClip.second, read1->entryQual.end());
				}
				mainQualPStore->resize(qualTmp->size());
				fastaPhredsToLog10Prob(qualTmp->size(), (unsigned char*)(&((*qualTmp)[0])), &((*mainQualPStore)[0]));
			//get the competing references
				std::pair<uintptr_t,uintptr_t> refAGot(std::max((intptr_t)0, breakInd-(intptr_t)readAlnSize), breakInd);
				std::pair<uintptr_t,uintptr_t> refBGot(breakInd, std::min(selRef->size(), breakInd+readAlnSize));
				refATmp->clear(); refATmp->insert(refATmp->end(), selRef->begin() + refAGot.first, selRef->begin() + refAGot.second);
				refBTmp->clear(); refBTmp->insert(refBTmp->end(), selRef->begin() + refBGot.first, selRef->begin() + refBGot.second);
			//likelihood of A
				PositionDependentCostKDTree* refACosts = &(rebaseCosts[threadInd]);
					refACosts->regionsRebased(selCost, refAGot.first, refAGot.second, -1, -1);
				if(selMang){
					PositionDependentQualityMangleSet* subMang = &(rebaseMangs[threadInd]);
						subMang->rebase(selMang, refAGot.first, refAGot.second);
					PositionDependentCostKDTree* tmpUse = &(mangleCosts[threadInd]);
						tmpUse->regionsQualityMangled(refACosts, subMang, qualTmp);
					refACosts = tmpUse;
				}
				refACosts->produceFromRegions();
				PositionDependentAffineGapLinearPairwiseAlignment* refAAln = &(workAlns[threadInd]);
				refAAln->changeProblem(2, refATmp, seqTmp, refACosts);
				refAAln->prepareAlignmentStructure();
				LinearPairwiseAlignmentIteration* refAIter = workIters[threadInd];
				if(!refAIter){
					workIters[threadInd] = refAAln->getIteratorToken();
					refAIter = workIters[threadInd];
				}
				int numScoreA = refAAln->findAlignmentScores(refAIter, mainScores->size(), &((*mainScores)[0]), worstScore, hotfuzz);
				packScoreSeen->resize(numScoreA);
				for(int i = 0; i<numScoreA; i++){ (*packScoreSeen)[i] = uptoCount; }
				refAAln->startFuzzyIteration(refAIter, (*mainScores)[numScoreA-1], hotfuzz, numScoreA);
				double refALPro = linearReferenceSourceProbabilityAffine(refAAln, &((*mainQualPStore)[0]), refAIter, lproGapOpen, lproGapExtend, numScoreA, &((*mainScores)[0]), uptoCount ? &((*packScoreSeen)[0]) : 0);
			//likelihood of B
				PositionDependentCostKDTree* refBCosts = &(rebaseCosts[threadInd]);
					refBCosts->regionsRebased(selCost, refBGot.first, refBGot.second, -1, -1);
				if(selMang){
					PositionDependentQualityMangleSet* subMang = &(rebaseMangs[threadInd]);
						subMang->rebase(selMang, refBGot.first, refBGot.second);
					PositionDependentCostKDTree* tmpUse = &(mangleCosts[threadInd]);
						tmpUse->regionsQualityMangled(refBCosts, subMang, qualTmp);
					refBCosts = tmpUse;
				}
				refBCosts->produceFromRegions();
				PositionDependentAffineGapLinearPairwiseAlignment* refBAln = &(workAlns[threadInd]);
				refBAln->changeProblem(2, refBTmp, seqTmp, refBCosts);
				refBAln->prepareAlignmentStructure();
				LinearPairwiseAlignmentIteration* refBIter = workIters[threadInd];
				if(!refBIter){
					workIters[threadInd] = refBAln->getIteratorToken();
					refBIter = workIters[threadInd];
				}
				int numScoreB = refBAln->findAlignmentScores(refBIter, mainScores->size(), &((*mainScores)[0]), worstScore, hotfuzz);
				packScoreSeen->resize(numScoreB);
				for(int i = 0; i<numScoreB; i++){ (*packScoreSeen)[i] = uptoCount; }
				refBAln->startFuzzyIteration(refBIter, (*mainScores)[numScoreB-1], hotfuzz, numScoreB);
				double refBLPro = linearReferenceSourceProbabilityAffine(refBAln, &((*mainQualPStore)[0]), refBIter, lproGapOpen, lproGapExtend, numScoreB, &((*mainScores)[0]), uptoCount ? &((*packScoreSeen)[0]) : 0);
			//make the decision
				canMoveRightOfLeft = (refALPro - refBLPro) > threshLR;
		}
	//final check
		switch(read1HPI - read2LPI){
			case 0:
				if(canMoveLeftOfRight || canMoveRightOfLeft){
					return 0;
				}
				break;
			case 1:
				if(canMoveLeftOfRight && canMoveRightOfLeft){
					return 0;
				}
				if(canMoveLeftOfRight && ((*selProbP)[read1HPI])){
					return 0;
				}
				if(canMoveRightOfLeft && ((*selProbP)[read2LPI])){
					return 0;
				}
				break;
			case 2:
				if(canMoveLeftOfRight && canMoveRightOfLeft){
					return 0; //already know it's problematic
				}
				break;
			default:
				std::cerr << "Da fuq?" << std::endl;
				return 0;
		}
		return 1;
}

ProsynarFilter* factoryProblematicRegionFilter(){
	return new ProblematicRegionFilter();
}
