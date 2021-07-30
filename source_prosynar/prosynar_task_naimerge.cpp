#include "prosynar_task_naimerge.h"

#include <math.h>
#include <string.h>
#include <algorithm>

#include "whodun_parse_seq.h"

void mergeIndelExpandedSequences(std::vector<int>* seqA, std::vector<double>* qualA, std::vector<int>* seqB, std::vector<double>* qualB, std::string* mergeSeq, std::vector<double>* mergeQual){
	for(uintptr_t i = 0; i<seqA->size(); i++){
		int sac = (*seqA)[i];
		double lqac = (*qualA)[i];
		int sbc = (*seqB)[i];
		double lqbc = (*qualB)[i];
		double errA = pow(10.0, lqac);
		double errB = pow(10.0, lqbc);
		if(lqac < lqbc){
			if(sac >= 0){
				mergeSeq->push_back(sac);
				if(sac != sbc){
					double numVal = errA * (1.0 - errB/3.0);
					double denVal = errA + errB - 4*errA*errB/3.0;
					mergeQual->push_back(log10(numVal/denVal));
				}
				else{
					double numVal = errA*errB / 3;
					double denVal = 1 - errA - errB + 4*errA*errB/3.0;
					mergeQual->push_back(log10(numVal/denVal));
				}
			}
		}
		else{
			if(sbc >= 0){
				mergeSeq->push_back(sbc);
				if(sac != sbc){
					double numVal = errB * (1.0 - errA/3.0);
					double denVal = errA + errB - 4*errA*errB/3.0;
					mergeQual->push_back(log10(numVal/denVal));
				}
				else{
					double numVal = errA*errB / 3;
					double denVal = 1 - errA - errB + 4*errA*errB/3.0;
					mergeQual->push_back(log10(numVal/denVal));
				}
			}
		}
	}
}

SimpleAlignMerger::SimpleAlignMerger(){
	softReclaim = false;
	reqOverlap = 1;
	qualmFile = 0;
	costReadFile = 0;
	myMainDoc = "prosynar -- Malign [OPTION]\nMerge by aligning the sequences.\nThe OPTIONS are:\n";
	myVersionDoc = "ProSynAr Malign 1.0";
	myCopyrightDoc = "Copyright (C) 2019 UNT HSC Center for Human Identification";
	ArgumentParserIntMeta overMeta("Overlap Threshold");
		addIntegerOption("--over", &reqOverlap, 0, "    Specify the number of bases of overlap to require.\n    --over 10\n", &overMeta);
	ArgumentParserBoolMeta promMeta("Reclaim Soft Clipped Bases");
		addBooleanFlag("--unclip", &softReclaim, 1, "    Reclaim bases that were originally clipped.\n", &promMeta);
	ArgumentParserStrMeta costMeta("Alignment Parameter File");
		costMeta.isFile = true;
		costMeta.fileExts.insert(".agc");
		addStringOption("--cost", &costReadFile, 0, "    Specify the alignment cost parameters.\n    --cost File.agc\n", &costMeta);
	ArgumentParserStrMeta qualmMeta("Biquality Mangle File");
		qualmMeta.isFile = true;
		qualmMeta.fileExts.insert(".bqualm");
		addStringOption("--bqualm", &qualmFile, 0, "    Specify how to modify alignment parameters using quality.\n    --bqualm File.bqualm\n", &qualmMeta);
}

SimpleAlignMerger::~SimpleAlignMerger(){
	for(uintptr_t i = 0; i<runIters.size(); i++){
		if(runIters[i]){ delete(runIters[i]); }
	}
	runIters.clear();
}

int SimpleAlignMerger::handleUnknownArgument(int argc, char** argv, std::ostream* helpOut){
	if(strcmp(argv[0],"--")==0){
		return 0;
	}
	argumentError.append("Unknown argument ");
	argumentError.append(argv[0]);
	return -1;
}

void SimpleAlignMerger::printExtraGUIInformation(std::ostream* toPrint){
	//nothing special
}

int SimpleAlignMerger::posteriorCheck(){
	if(reqOverlap < 0){
		argumentError = "Overlap threshold must be non-negative.";
		return 1;
	}
	if(costReadFile == 0){
		argumentError = "Alignment parameter file required.";
		return 1;
	}
	return 0;
}

void SimpleAlignMerger::initialize(ProsynarArgumentParser* baseArgs){
	saveArgs = baseArgs;
	std::string fileConts;
	FileInStream redFile(costReadFile);
	readStream(&redFile, &fileConts);
	mergeCost.parseSpecification(fileConts.c_str(), fileConts.c_str() + fileConts.size());
	bigCost.regionsUniform(&mergeCost);
	bigCost.produceFromRegions();
	if(qualmFile){
		std::string qfileConts;
		FileInStream redQFile(qualmFile);
		readStream(&redQFile, &qfileConts);
		bigMang.parseQualityMangleSet(qfileConts.c_str(), qfileConts.c_str() + qfileConts.size());
	}
	cigLocSet.resize(baseArgs->numThread);
	costSet.resize(baseArgs->numThread);
	seqASet.resize(baseArgs->numThread);
	seqqASet.resize(baseArgs->numThread);
	seqqdASet.resize(baseArgs->numThread);
	seqBSet.resize(baseArgs->numThread);
	seqqBSet.resize(baseArgs->numThread);
	seqqdBSet.resize(baseArgs->numThread);
	seqILSet.resize(baseArgs->numThread);
	seqdILSet.resize(baseArgs->numThread);
	seqIRSet.resize(baseArgs->numThread);
	seqdIRSet.resize(baseArgs->numThread);
	saveAlns.resize(baseArgs->numThread);
	LinearPairwiseAlignmentIteration* nullIt = 0;
	runIters.insert(runIters.end(), baseArgs->numThread, nullIt);
}

int SimpleAlignMerger::mergePair(int threadInd, CRBSAMFileContents* read1, CRBSAMFileContents* read2, std::string* mergeSeq, std::vector<double>* mergeQual, std::string* errRep){
	const char* read1SStart = &(read1->entrySeq[0]);
	const char* read1QStart = &(read1->entryQual[0]);
	uintptr_t read1SLen = read1->entrySeq.size();
	const char* read2SStart = &(read2->entrySeq[0]);
	const char* read2QStart = &(read2->entryQual[0]);
	uintptr_t read2SLen = read2->entrySeq.size();
	if(!softReclaim){
		std::vector<intptr_t>* cigVec = &(cigLocSet[threadInd]);
		//expand the cigars (really want the soft clips)
		std::pair<uintptr_t,uintptr_t> read1SClip;
		try{
			cigVec->clear();
			read1SClip = cigarStringToReferencePositions(read1->entryPos, &(read1->entryCigar), cigVec);
		}catch(std::exception& err){
			errRep->insert(errRep->end(),read1->entryName.begin(), read1->entryName.end());
			errRep->push_back(':');
			errRep->push_back(' ');
			errRep->append(err.what());
			return -1;
		}
		std::pair<uintptr_t,uintptr_t> read2SClip;
		try{
			cigVec->clear();
			read2SClip = cigarStringToReferencePositions(read2->entryPos, &(read2->entryCigar), cigVec);
		}catch(std::exception& err){
			errRep->insert(errRep->end(),read2->entryName.begin(), read2->entryName.end());
			errRep->push_back(':');
			errRep->push_back(' ');
			errRep->append(err.what());
			return -1;
		}
		//do the clipping
		read1SStart += read1SClip.first;
		read1QStart += read1SClip.first;
		read1SLen -= (read1SClip.first + read1SClip.second);
		read2SStart += read2SClip.first;
		read2QStart += read2SClip.first;
		read2SLen -= (read2SClip.first + read2SClip.second);
	}
	//quick abandon if too short
		if((read1SLen < (uintptr_t)reqOverlap) || (read2SLen < (uintptr_t)reqOverlap)){
			return 1;
		}
	//get the stuff
		std::string* seqA = &(seqASet[threadInd]);
		std::vector<char>* seqAQ = &(seqqASet[threadInd]);
		std::vector<double>* seqAQD = &(seqqdASet[threadInd]);
			seqA->clear(); seqA->insert(seqA->end(), read1SStart, read1SStart + read1SLen);
			seqAQ->clear(); seqAQ->insert(seqAQ->end(), read1QStart, read1QStart + read1SLen);
			seqAQD->resize(read1SLen); fastaPhredsToLog10Prob(read1SLen, (const unsigned char*)read1QStart, &((*seqAQD)[0]));
		std::string* seqB = &(seqBSet[threadInd]);
		std::vector<char>* seqBQ = &(seqqBSet[threadInd]);
		std::vector<double>* seqBQD = &(seqqdBSet[threadInd]);
			seqB->clear(); seqB->insert(seqB->end(), read2SStart, read2SStart + read2SLen);
			seqBQ->clear(); seqBQ->insert(seqBQ->end(), read2QStart, read2QStart + read2SLen);
			seqBQD->resize(read2SLen); fastaPhredsToLog10Prob(read2SLen, (const unsigned char*)read2QStart, &((*seqBQD)[0]));
	//mangle the alignment parameters
		PositionDependentCostKDTree* useCost = &bigCost;
		if(qualmFile){
			useCost = &(costSet[threadInd]);
			useCost->regionsBiQualityMangled(&bigCost, &bigMang, seqAQ, seqBQ);
			useCost->produceFromRegions();
		}
	//do an alignment
		PositionDependentAffineGapLinearPairwiseAlignment* curAln = &(saveAlns[threadInd]);
		curAln->changeProblem(2, seqA, seqB, useCost);
		curAln->prepareAlignmentStructure();
		LinearPairwiseAlignmentIteration* curIter = runIters[threadInd];
		if(!curIter){
			runIters[threadInd] = curAln->getIteratorToken();
			curIter = runIters[threadInd];
		}
		curAln->startOptimalIteration(curIter);
		int wasAln = curIter->getNextAlignment();
		while(wasAln){
			if(curIter->aInds.size() > (uintptr_t)reqOverlap){
				break;
			}
			wasAln = curIter->getNextAlignment();
		}
		if(!wasAln){
			return 1;
		}
		uintptr_t alnLen = curIter->aInds.size();
	//expand out indels
		std::vector<int>* seqIL = &(seqILSet[threadInd]);
		std::vector<double>* seqdIL = &(seqdILSet[threadInd]);
		seqIL->clear(); seqdIL->clear();
		std::vector<int>* seqIR = &(seqIRSet[threadInd]);
		std::vector<double>* seqdIR = &(seqdIRSet[threadInd]);
		seqIR->clear(); seqdIR->clear();
		for(uintptr_t i = 1; i<alnLen; i++){
			if(curIter->aInds[i] != curIter->aInds[i-1]){
				if(curIter->bInds[i] != curIter->bInds[i-1]){
					//consume both
					seqIL->push_back((*seqA)[curIter->aInds[i-1]]);
					seqdIL->push_back((*seqAQD)[curIter->aInds[i-1]]);
					seqIR->push_back((*seqB)[curIter->bInds[i-1]]);
					seqdIR->push_back((*seqBQD)[curIter->bInds[i-1]]);
				}
				else{
					//insert a, skip b
					seqIL->push_back((*seqA)[curIter->aInds[i-1]]);
					seqdIL->push_back((*seqAQD)[curIter->aInds[i-1]]);
					seqIR->push_back(-1);
					seqdIR->push_back(0.0);
				}
			}
			else{
				//insert b, skip a
				seqIL->push_back(-1);
				seqdIL->push_back(0.0);
				seqIR->push_back((*seqB)[curIter->bInds[i-1]]);
				seqdIR->push_back((*seqBQD)[curIter->bInds[i-1]]);
			}
		}
	//fill in any indel qualities
		uintptr_t ai = 0;
		while(ai < seqIL->size()){
			if((*seqIL)[ai] >= 0){
				ai++; continue;
			}
			uintptr_t nai = ai + 1;
			while((nai < seqIL->size()) && ((*seqIL)[nai] < 0)){
				nai++;
			}
			double winQual;
			if(ai && (nai < seqIL->size())){
				double leftQ = (*seqdIL)[ai-1];
				double rightQ = (*seqdIL)[nai];
				winQual = log10((pow(10.0, leftQ) + pow(10.0, rightQ))/2.0);
			}
			else if(ai){ winQual = (*seqdIL)[ai-1]; }
			else if(nai < seqIL->size()){ winQual = (*seqdIL)[nai]; }
			else{ winQual = -1.0 / 0.0; }
			for(uintptr_t ci = ai; ci < nai; ci++){
				(*seqdIL)[ci] = winQual;
			}
			ai = nai;
		}
		uintptr_t bi = 0;
		while(bi < seqIL->size()){
			if((*seqIL)[bi] >= 0){
				bi++; continue;
			}
			uintptr_t nbi = bi + 1;
			while((nbi < seqIL->size()) && ((*seqIL)[nbi] < 0)){
				nbi++;
			}
			double winQual;
			if(bi && (nbi < seqIL->size())){
				double leftQ = (*seqdIL)[bi-1];
				double rightQ = (*seqdIL)[nbi];
				winQual = log10((pow(10.0, leftQ) + pow(10.0, rightQ))/2.0);
			}
			else if(bi){ winQual = (*seqdIL)[bi-1]; }
			else if(nbi < seqIL->size()){ winQual = (*seqdIL)[nbi]; }
			else{ winQual = -1.0 / 0.0; }
			for(uintptr_t ci = bi; ci < nbi; ci++){
				(*seqdIL)[ci] = winQual;
			}
			bi = nbi;
		}
	//do the merge
		mergeSeq->clear(); mergeQual->clear();
		//prefix
		if(curIter->aInds[0] != 0){
			mergeSeq->insert(mergeSeq->end(), seqA->begin(), seqA->begin() + curIter->aInds[0]);
			mergeQual->insert(mergeQual->end(), seqAQD->begin(), seqAQD->begin() + curIter->aInds[0]);
		}
		else{
			mergeSeq->insert(mergeSeq->end(), seqB->begin(), seqB->begin() + curIter->bInds[0]);
			mergeQual->insert(mergeQual->end(), seqBQD->begin(), seqBQD->begin() + curIter->bInds[0]);
		}
		//overlap
		mergeIndelExpandedSequences(seqIL, seqdIL, seqIR, seqdIR, mergeSeq, mergeQual);
		//suffix
		if((uintptr_t)(curIter->aInds[alnLen-1]) < seqA->size()){
			mergeSeq->insert(mergeSeq->end(), seqA->begin() + curIter->aInds[alnLen-1], seqA->end());
			mergeQual->insert(mergeQual->end(), seqAQD->begin() + curIter->aInds[alnLen-1], seqAQD->end());
		}
		else{
			mergeSeq->insert(mergeSeq->end(), seqB->begin() + curIter->bInds[alnLen-1], seqB->end());
			mergeQual->insert(mergeQual->end(), seqBQD->begin() + curIter->bInds[alnLen-1], seqBQD->end());
		}
	return 0;
}

ProsynarMerger* factorySimpleAlignMerger(){
	return new SimpleAlignMerger();
}

