#include "prosynar_task_pear.h"

#include <math.h>
#include <float.h>
#include <string.h>
#include <iostream>
#include <algorithm>

#include "whodun_probutil.h"
#include "whodun_parse_seq.h"
#include "prosynar_task_naimerge.h"

PEARMerger::PEARMerger(){
	softReclaim = false;
	reqOverlap = 1;
	sigLevel = 0.01;
	origRecipe = false;
	matchPoints = 1.0;
	mismatchPoints = -1.0;
	testObsOver = false;
	myMainDoc = "prosynar -- Mpear [OPTION]\nMerge by sliding the sequences past each other\nlooking for the highest score.\nThe OPTIONS are:\n";
	myVersionDoc = "ProSynAr Mpear 1.0";
	myCopyrightDoc = "Copyright (C) 2020 UNT HSC Center for Human Identification\nCITE: Zhang, J., Kobert, K., Flouri, T., & Stamatakis, A. (2014). PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics, 30(5), 614-620.";
	ArgumentParserIntMeta overMeta("Overlap Threshold");
		addIntegerOption("--over", &reqOverlap, 0, "    Specify the number of bases of overlap to require.\n    --over 1\n", &overMeta);
	ArgumentParserBoolMeta promMeta("Reclaim Soft Clipped Bases");
		addBooleanFlag("--unclip", &softReclaim, 1, "    Reclaim bases that were originally clipped.\n", &promMeta);
	ArgumentParserBoolMeta qualMeta("Original Quality Method");
		addBooleanFlag("--pearqual", &origRecipe, 1, "    Use the method specified in the original PEAR paper.\n    Defaults to the method by Edgar and Flyvbjerg (2015).\n", &qualMeta);
	ArgumentParserFltMeta sigMeta("Significance Level");
		addFloatOption("--alpha", &sigLevel, 0, "    Set the significance level of the score test.\n    --alpha 0.01\n", &sigMeta);
	ArgumentParserFltMeta matMeta("Match Score");
		addFloatOption("--ms", &matchPoints, 0, "    Points for matching bases.\n    --ms 1.0\n", &matMeta);
	ArgumentParserFltMeta mmatMeta("Mismatch Score");
		addFloatOption("--xs", &mismatchPoints, 0, "    Points for non-matching bases.\n    --xs -1.0\n", &mmatMeta);
	ArgumentParserBoolMeta mapMeta("Maximal Accepted Probability");
		addBooleanFlag("--map", &testObsOver, 1, "    Do the significance test against the observed overlap.\n", &mapMeta);
}

PEARMerger::~PEARMerger(){
}

int PEARMerger::handleUnknownArgument(int argc, char** argv, std::ostream* helpOut){
	if(strcmp(argv[0],"--")==0){
		return 0;
	}
	argumentError.append("Unknown argument ");
	argumentError.append(argv[0]);
	return -1;
}

void PEARMerger::printExtraGUIInformation(std::ostream* toPrint){
	//nothing special
}

int PEARMerger::posteriorCheck(){
	if(reqOverlap <= 0){
		argumentError = "Overlap threshold must be positive.";
		return 1;
	}
	if(sigLevel < 0){
		argumentError = "Significance level must be non-negative.";
		return 1;
	}
	if(sigLevel > 1){
		argumentError = "Significance level must be less than 1.";
		return 1;
	}
	return 0;
}

void PEARMerger::initialize(ProsynarArgumentParser* baseArgs){
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
 * @param matPts The points for a match.
 * @param mmatPts The points for a mismatch.
 * @return The number of bases of overlap, and the score.
 */
std::pair<uintptr_t,double> pearFindBestOverlap(std::vector<char>* seqL, std::vector<double>* qualL, std::vector<char>* seqR, std::vector<double>* qualR, double matPts, double mmatPts){
	double winScore = -1.0/0.0;
	uintptr_t winOver = 0;
	uintptr_t maxPos = std::min(seqL->size(), seqR->size());
	//count the base frequencies
	double allFreqs[256];
	memset(allFreqs, 0, 256*sizeof(double));
	for(uintptr_t i = 0; i<seqL->size(); i++){
		allFreqs[(*seqL)[i] & 0x00FF] += 1.0;
	}
	for(uintptr_t i = 0; i<seqR->size(); i++){
		allFreqs[(*seqR)[i] & 0x00FF] += 1.0;
	}
	uintptr_t totNBase = seqL->size() + seqR->size();
	for(uintptr_t i = 0; i<256; i++){
		allFreqs[i] = allFreqs[i] / totNBase;
	}
	//get the best score
	for(uintptr_t col = 1; col <= maxPos; col++){
		uintptr_t sli0 = seqL->size() - col;
		double curScore = 0.0;
		for(uintptr_t i = 0; i<col; i++){
			unsigned char baseL = 0x00FF & ((*seqL)[sli0 + i]);
			unsigned char baseR = 0x00FF & ((*seqR)[i]);
			double errL = pow(10.0, (*qualL)[sli0 + i]);
			double errR = pow(10.0, (*qualR)[i]);
			if(baseL == baseR){
				double curNum = 0.0; double curDen = 0.0;
				for(uintptr_t i = 0; i<256; i++){
					if(i == baseL){ continue; }
					double curBPro = allFreqs[i];
					curNum += (curBPro*curBPro);
					curDen += curBPro;
				}
				curDen = curDen*curDen;
				curScore += matPts*((1-errL)*(1-errR) + errL*errR*curNum/curDen);
			}
			else{
				double curNumA = allFreqs[baseR]; double curDenA = 0.0;
				double curNumB = allFreqs[baseL]; double curDenB = 0.0;
				double curNumC = 0.0; double curDenC = 0.0;
				for(uintptr_t i = 0; i<256; i++){
					if(i != baseL){ curDenA += allFreqs[i]; }
					if(i != baseR){ curDenB += allFreqs[i]; }
					if((i!=baseL) && (i!=baseR)){
						double curBPro = allFreqs[i];
						curNumC += (curBPro*curBPro);
						curDenC += curBPro;
					}
				}
				curDenC = curDenC*curDenC;
				double curPro = (1-errR)*errL*curNumA/curDenA;
					curPro += (1-errL)*errR*curNumB/curDenB;
					curPro += errL*errR*curNumC/curDenC;
				curScore += mmatPts*(1.0-curPro);
			}
		}
		if(curScore > winScore){
			winScore = curScore;
			winOver = col;
		}
	}
	return std::pair<uintptr_t,double>(winOver,winScore);
}

int PEARMerger::mergePair(int threadInd, CRBSAMFileContents* read1, CRBSAMFileContents* read2, std::string* mergeSeq, std::vector<double>* mergeQual, std::string* errRep){
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
		std::pair<uintptr_t,double> resAB = pearFindBestOverlap(seqA, seqAQD, seqB, seqBQD, matchPoints, mismatchPoints);
		std::pair<uintptr_t,double> resBA = pearFindBestOverlap(seqB, seqBQD, seqA, seqAQD, matchPoints, mismatchPoints);
		std::pair<uintptr_t,double> resAR = pearFindBestOverlap(seqA, seqAQD, seqBR, seqBRQD, matchPoints, mismatchPoints);
		std::pair<uintptr_t,double> resRA = pearFindBestOverlap(seqBR, seqBRQD, seqA, seqAQD, matchPoints, mismatchPoints);
	//winner winner
		std::vector<char>* winL = seqA; /*std::vector<char>* winLQ = seqAQ;*/ std::vector<double>* winLQD = seqAQD;
		std::vector<char>* winR = seqB; /*std::vector<char>* winRQ = seqBQ;*/ std::vector<double>* winRQD = seqBQD;
		std::pair<uintptr_t,double> winRes = resAB;
		if(resBA.second > winRes.second){
			winL = seqB; /*winLQ = seqBQ;*/ winLQD = seqBQD;
			winR = seqA; /*winRQ = seqAQ;*/ winRQD = seqAQD;
			winRes = resBA;
		}
		if(resAR.second > winRes.second){
			winL = seqA; /*winLQ = seqAQ;*/ winLQD = seqAQD;
			winR = seqBR; /*winRQ = seqBRQ;*/ winRQD = seqBRQD;
			winRes = resAR;
		}
		if(resRA.second > winRes.second){
			winL = seqBR; /*winLQ = seqBRQ;*/ winLQD = seqBRQD;
			winR = seqA; /*winRQ = seqAQ;*/ winRQD = seqAQD;
			winRes = resRA;
		}
	//if too short, stop
		if(winRes.first < (uintptr_t)reqOverlap){
			return 1;
		}
	//count the base frequencies
		double allFreqs[256];
		memset(allFreqs, 0, 256*sizeof(double));
		for(uintptr_t i = 0; i<winL->size(); i++){
			allFreqs[(*winL)[i] & 0x00FF] += 1.0;
		}
		for(uintptr_t i = 0; i<winR->size(); i++){
			allFreqs[(*winR)[i] & 0x00FF] += 1.0;
		}
		uintptr_t totNBase = winL->size() + winR->size();
		for(uintptr_t i = 0; i<256; i++){
			allFreqs[i] = allFreqs[i] / totNBase;
		}
	//score p-test
		double ranMatP = 0.0;
		for(uintptr_t i = 0; i<256; i++){
			ranMatP += (allFreqs[i]*allFreqs[i]);
		}
		double ranMMP = 1.0 - ranMatP;
		ranMatP = log(ranMatP);
		ranMMP = log(ranMMP);
		double lnOf10 = log(10.0);
		uintptr_t curOver = testObsOver ? winRes.first : reqOverlap;
		uintptr_t maxOver = 2 * std::max(winL->size(), winR->size());
		double totL10Pro = 0.0;
		while(curOver < maxOver){
			double addToD = ceil((winRes.second - mismatchPoints*curOver) / (matchPoints - mismatchPoints)) - 1;
			if(addToD < 0){ curOver++; continue; }
			uintptr_t addTo = (uintptr_t)addToD;
			ProbabilitySummation curSum;
			for(uintptr_t k = 0; k<=addTo; k++){
				double curLnEnt = logGamma(curOver+1);
				curLnEnt -= logGamma(k + 1);
				curLnEnt -= logGamma(curOver - k + 1);
				curLnEnt += k*ranMatP;
				curLnEnt += (curOver - k)*ranMMP;
				curSum.addNextLogProb(curLnEnt / lnOf10);
			}
			double curL10Mul = curSum.getFinalLogSum();
			/* Need a way to quit early, but this doesn't work.
			if(fabs(curL10Mul) < DBL_EPSILON){
				break;
			}
			*/
			totL10Pro += curL10Mul;
			curOver++;
		}
		double pvalue = 1.0 - pow(10.0, totL10Pro);
	//decide
		if(pvalue > sigLevel){
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
					mergeQual->push_back( (optL == optR) ? (optLQ + optRQ) : optLQ);
				}
				else{
					mergeSeq->push_back(optR);
					mergeQual->push_back( (optL == optR) ? (optLQ + optRQ) : optRQ);
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

ProsynarMerger* factoryPearMerger(){
	return new PEARMerger();
}


