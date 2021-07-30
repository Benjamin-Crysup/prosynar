#include "prosynar_task_flash.h"

#include <math.h>
#include <string.h>
#include <algorithm>

#include "whodun_parse_seq.h"
#include "prosynar_task_naimerge.h"

FLASHMerger::FLASHMerger(){
	softReclaim = false;
	reqOverlap = 10;
	maxExpOverlap = 70;
	origRecipe = false;
	worstMR = 0.1;
	myMainDoc = "prosynar -- Mflash [OPTION]\nMerge by sliding the sequences past each other.\nThe OPTIONS are:\n";
	myVersionDoc = "ProSynAr Mflash 1.0";
	myCopyrightDoc = "Copyright (C) 2019 UNT HSC Center for Human Identification\nCITE: Magoc, T., & Salzberg, S. L. (2011). FLASH: fast length adjustment of short reads to improve genome assemblies. Bioinformatics, 27(21), 2957-2963.";
	ArgumentParserIntMeta overMeta("Overlap Threshold");
		addIntegerOption("--over", &reqOverlap, 0, "    Specify the number of bases of overlap to require.\n    --over 10\n", &overMeta);
	ArgumentParserBoolMeta promMeta("Reclaim Soft Clipped Bases");
		addBooleanFlag("--unclip", &softReclaim, 1, "    Reclaim bases that were originally clipped.\n", &promMeta);
	ArgumentParserIntMeta maxoverMeta("Max Overlap Expected");
		addIntegerOption("--maxexp", &maxExpOverlap, 0, "    Specify the expected number of bases of overlap.\n    --maxexp 70\n", &maxoverMeta);
	ArgumentParserBoolMeta qualMeta("Original Quality Method");
		addBooleanFlag("--flashqual", &origRecipe, 1, "    Use the method specified in the original FLASH paper.\n    Defaults to the method by Edgar and Flyvbjerg (2015).\n", &qualMeta);
	ArgumentParserFltMeta misMeta("Mismatch Threshold");
		addFloatOption("--misrat", &worstMR, 0, "    Specify the worst mismatch ratio to accept.\n    --misrat 1.0\n", &misMeta);
}

FLASHMerger::~FLASHMerger(){
}

int FLASHMerger::handleUnknownArgument(int argc, char** argv, std::ostream* helpOut){
	if(strcmp(argv[0],"--")==0){
		return 0;
	}
	argumentError.append("Unknown argument ");
	argumentError.append(argv[0]);
	return -1;
}

void FLASHMerger::printExtraGUIInformation(std::ostream* toPrint){
	//nothing special
}

int FLASHMerger::posteriorCheck(){
	if(reqOverlap <= 0){
		argumentError = "Overlap threshold must be positive.";
		return 1;
	}
	if(maxExpOverlap <= 0){
		argumentError = "Expected overlap must be positive.";
		return 1;
	}
	if(worstMR < 0){
		argumentError = "Mismatch ratio threshold must be non-negative.";
		return 1;
	}
	if(worstMR > 1){
		argumentError = "Mismatch ratio threshold must be less than 1.";
		return 1;
	}
	return 0;
}

void FLASHMerger::initialize(ProsynarArgumentParser* baseArgs){
	saveArgs = baseArgs;
	cigLocSet.resize(baseArgs->numThread);
	seqASet.resize(baseArgs->numThread);
	seqqASet.resize(baseArgs->numThread);
	seqqdASet.resize(baseArgs->numThread);
	seqBSet.resize(baseArgs->numThread);
	seqqBSet.resize(baseArgs->numThread);
	seqqdBSet.resize(baseArgs->numThread);
	seqBRSet.resize(baseArgs->numThread);
	seqqBRSet.resize(baseArgs->numThread);
	seqqdBRSet.resize(baseArgs->numThread);
	seqILSet.resize(baseArgs->numThread);
	seqdILSet.resize(baseArgs->numThread);
	seqIRSet.resize(baseArgs->numThread);
	seqdIRSet.resize(baseArgs->numThread);
}

/**
 * Find the best overlap using the flash method.
 * @param seqL The left sequence.
 * @param qualL The left quality;
 * @param seqR The right sequence.
 * @param qualR The right quality.
 * @param minOver The minimum possible overlap.
 * @param maxExp THe maximum expected overlap.
 */
std::pair<uintptr_t,double> flashFindBestOverlap(std::vector<char>* seqL, std::vector<char>* qualL, std::vector<char>* seqR, std::vector<char>* qualR, uintptr_t minOver, uintptr_t maxExp){
	double winRat = 1.0/0.0;
	uintptr_t winOver = 0;
	double winAvgQ = 0.0;
	uintptr_t maxPos = std::min(seqL->size(), seqR->size());
	for(uintptr_t col = minOver; col <= maxPos; col++){
		uintptr_t sli0 = seqL->size() - col;
		uintptr_t numMiss = 0;
		double curAvgQ = 0.0;
		for(uintptr_t i = 0; i<col; i++){
			if((*seqL)[sli0 + i] != (*seqR)[i]){
				curAvgQ += (*qualL)[sli0+i] + (*qualR)[i];
				numMiss++;
			}
		}
		double curRat = numMiss;
		if(col > maxExp){
			curRat = curRat / maxExp;
		}
		else{
			curRat = curRat / col;
		}
		curAvgQ = curAvgQ / (2*numMiss);
		if((curRat < winRat) || ((curRat == winRat) && (curAvgQ < winAvgQ))){
			winRat = curRat;
			winAvgQ = curAvgQ;
			winOver = col;
		}
	}
	return std::pair<uintptr_t,double>(winOver,winRat);
}

int FLASHMerger::mergePair(int threadInd, CRBSAMFileContents* read1, CRBSAMFileContents* read2, std::string* mergeSeq, std::vector<double>* mergeQual, std::string* errRep){
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
		std::vector<char>* seqA = &(seqASet[threadInd]);
		std::vector<char>* seqAQ = &(seqqASet[threadInd]);
		std::vector<double>* seqAQD = &(seqqdASet[threadInd]);
			seqA->clear(); seqA->insert(seqA->end(), read1SStart, read1SStart + read1SLen);
			seqAQ->clear(); seqAQ->insert(seqAQ->end(), read1QStart, read1QStart + read1SLen);
			seqAQD->resize(read1SLen); fastaPhredsToLog10Prob(read1SLen, (const unsigned char*)read1QStart, &((*seqAQD)[0]));
		std::vector<char>* seqB = &(seqBSet[threadInd]);
		std::vector<char>* seqBQ = &(seqqBSet[threadInd]);
		std::vector<double>* seqBQD = &(seqqdBSet[threadInd]);
			seqB->clear(); seqB->insert(seqB->end(), read2SStart, read2SStart + read2SLen);
			seqBQ->clear(); seqBQ->insert(seqBQ->end(), read2QStart, read2QStart + read2SLen);
			seqBQD->resize(read2SLen); fastaPhredsToLog10Prob(read2SLen, (const unsigned char*)read2QStart, &((*seqBQD)[0]));
		std::vector<char>* seqBR = &(seqBRSet[threadInd]);
		std::vector<char>* seqBRQ = &(seqqBRSet[threadInd]);
		std::vector<double>* seqBRQD = &(seqqdBRSet[threadInd]);
			seqBR->clear(); seqBR->insert(seqBR->end(), read2SStart, read2SStart + read2SLen);
			seqBRQ->clear(); seqBRQ->insert(seqBRQ->end(), read2QStart, read2QStart + read2SLen);
			seqBRQD->clear(); seqBRQD->insert(seqBRQD->end(), seqBQD->begin(), seqBQD->end());
			sequenceReverseCompliment(read2SLen, &((*seqBR)[0]), &((*seqBRQD)[0]));
			std::reverse(seqBRQ->begin(), seqBRQ->end());
	//start walking
		std::pair<uintptr_t,double> resAB = flashFindBestOverlap(seqA, seqAQ, seqB, seqBQ, reqOverlap, maxExpOverlap);
		std::pair<uintptr_t,double> resBA = flashFindBestOverlap(seqB, seqBQ, seqA, seqAQ, reqOverlap, maxExpOverlap);
		std::pair<uintptr_t,double> resAR = flashFindBestOverlap(seqA, seqAQ, seqBR, seqBRQ, reqOverlap, maxExpOverlap);
		std::pair<uintptr_t,double> resRA = flashFindBestOverlap(seqBR, seqBRQ, seqA, seqAQ, reqOverlap, maxExpOverlap);
	//winner winner
		std::vector<char>* winL = seqA; /*std::vector<char>* winLQ = seqAQ;*/ std::vector<double>* winLQD = seqAQD;
		std::vector<char>* winR = seqB; /*std::vector<char>* winRQ = seqBQ;*/ std::vector<double>* winRQD = seqBQD;
		std::pair<uintptr_t,double> winRes = resAB;
		if(resBA.second < winRes.second){
			winL = seqB; /*winLQ = seqBQ;*/ winLQD = seqBQD;
			winR = seqA; /*winRQ = seqAQ;*/ winRQD = seqAQD;
			winRes = resBA;
		}
		if(resAR.second < winRes.second){
			winL = seqA; /*winLQ = seqAQ;*/ winLQD = seqAQD;
			winR = seqBR; /*winRQ = seqBRQ;*/ winRQD = seqBRQD;
			winRes = resAR;
		}
		if(resRA.second < winRes.second){
			winL = seqBR; /*winLQ = seqBRQ;*/ winLQD = seqBRQD;
			winR = seqA; /*winRQ = seqAQ;*/ winRQD = seqAQD;
			winRes = resRA;
		}
		if(winRes.second > worstMR){
			return 1;
		}
	//chicken dinner
		mergeSeq->clear(); mergeQual->clear();
		if(origRecipe){
			uintptr_t sli0 = winL->size() - winRes.first;
			for(uintptr_t i = 0; i<sli0; i++){
				mergeSeq->push_back((*winL)[i]);
				mergeQual->push_back((*winLQD)[i]);
			}
			for(uintptr_t i = 0; i<winRes.first; i++){
				char optL = (*winL)[sli0+i];
				double optLQ = (*winLQD)[sli0+i];
				char optR = (*winR)[i];
				double optRQ = (*winRQD)[i];
				if(optLQ < optRQ){
					mergeSeq->push_back(optL);
					mergeQual->push_back( (optL == optR) ? optLQ : -0.2);
				}
				else{
					mergeSeq->push_back(optR);
					mergeQual->push_back( (optL == optR) ? optRQ : -0.2);
				}
			}
			for(uintptr_t i = winRes.first; i<winR->size(); i++){
				mergeSeq->push_back((*winR)[i]);
				mergeQual->push_back((*winRQD)[i]);
			}
		}
		else{
			uintptr_t sli0 = winL->size() - winRes.first;
			mergeSeq->insert(mergeSeq->end(), winL->begin(), winL->begin() + sli0);
			mergeQual->insert(mergeQual->end(), winLQD->begin(), winLQD->begin() + sli0);
			std::vector<int>* seqIL = &(seqILSet[threadInd]);
			std::vector<double>* seqdIL = &(seqdILSet[threadInd]);
			seqIL->clear(); seqdIL->clear();
			std::vector<int>* seqIR = &(seqIRSet[threadInd]);
			std::vector<double>* seqdIR = &(seqdIRSet[threadInd]);
			seqIR->clear(); seqdIR->clear();
			for(uintptr_t i = 0; i<winRes.first; i++){
				seqIL->push_back(0x00FF & (*winL)[sli0+i]);
				seqdIL->push_back((*winLQD)[sli0+i]);
				seqIR->push_back(0x00FF & (*winR)[i]);
				seqdIR->push_back((*winRQD)[i]);
			}
			mergeIndelExpandedSequences(seqIL, seqdIL, seqIR, seqdIR, mergeSeq, mergeQual);
			mergeSeq->insert(mergeSeq->end(), winR->begin() + winRes.first, winR->begin() + winR->size());
			mergeQual->insert(mergeQual->end(), winRQD->begin() + winRes.first, winRQD->begin() + winRQD->size());
		}
	return 0;
}

ProsynarMerger* factoryFlashMerger(){
	return new FLASHMerger();
}

