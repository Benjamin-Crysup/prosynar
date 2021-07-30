
#include <string.h>

#include "whodun_thread.h"
#include "whodun_oshook.h"
#include "whodun_parse_seq.h"
#include "whodun_genome_paired.h"

#include "prosynar_task.h"

/**A pair to try to merge.*/
class MergeAttemptTask{
public:
	/**The main alignment entry.*/
	CRBSAMFileContents* mainEnt;
	/**The pair entry.*/
	CRBSAMFileContents* pairEnt;
};

/**A merged sequence to output.*/
class MergeSequenceData{
public:
	/**The name of the sequence.*/
	std::string seqName;
	/**The merged sequence.*/
	std::string seqSeq;
	/**The merged qualities.*/
	std::vector<double> seqQuals;
	/**The original main alignment entry.*/
	CRBSAMFileContents* mainEnt;
	/**The original pair entry.*/
	CRBSAMFileContents* pairEnt;
};

/**Merge sequence thread arguments.*/
typedef struct{
	/**The index of this thread.*/
	int threadInd;
	/**Arguments to ProSynAr.*/
	ProsynarArgumentParser* argsP;
	/**Get tasks.*/
	ThreadProdComCollector<MergeAttemptTask>* getPCC;
	/**Return entry data.*/
	ThreadsafeReusableContainerCache<CRBSAMFileContents>* entC;
	/**Output failed.*/
	ThreadProdComCollector<CRBSAMFileContents>* failPCC;
	/**Output results to write.*/
	ThreadProdComCollector<MergeSequenceData>* goodPCC;
} MergeThreadArgs;
/**
 * Try to merge sequences.
 * @param tmpArg The arguments.
 */
void attemptMerging(void* tmpArg){
	MergeThreadArgs* myArgs = (MergeThreadArgs*)tmpArg;
	int myInd = myArgs->threadInd;
	std::string tmpErr;
	ProsynarArgumentParser* argsP = myArgs->argsP;
	MergeAttemptTask* anyRes = myArgs->getPCC->getThing();
	while(anyRes){
		tmpErr.clear();
		int mergeGreen = 1;
		for(uintptr_t i = 0; i<argsP->useFilters.size(); i++){
			int filtR = argsP->useFilters[i]->filterPair(myInd, anyRes->mainEnt, anyRes->pairEnt, &tmpErr);
			if(filtR < 0){
				mergeGreen = 0;
				lockMutex(argsP->errLock);
					std::cerr << tmpErr << std::endl;
				unlockMutex(argsP->errLock);
				break;
			}
			else if(filtR == 0){
				mergeGreen = 0;
				break;
			}
		}
		if(mergeGreen){
			MergeSequenceData* mrgInt = myArgs->goodPCC->taskCache.alloc();
			mrgInt->seqName.clear();
			mrgInt->seqSeq.clear();
			mrgInt->seqQuals.clear();
			mrgInt->seqName.insert(mrgInt->seqName.end(), anyRes->mainEnt->entryName.begin(), anyRes->mainEnt->entryName.end());
			mrgInt->mainEnt = anyRes->mainEnt;
			mrgInt->pairEnt = anyRes->pairEnt;
			int mergR = argsP->useMerger->mergePair(myInd, anyRes->mainEnt, anyRes->pairEnt, &(mrgInt->seqSeq), &(mrgInt->seqQuals), &tmpErr);
			if(mergR){
				mergeGreen = 0;
				if(mergR < 0){
					lockMutex(argsP->errLock);
						std::cerr << tmpErr << std::endl;
					unlockMutex(argsP->errLock);
				}
				myArgs->goodPCC->taskCache.dealloc(mrgInt);
			}
			else{
				myArgs->goodPCC->addThing(mrgInt);
			}
		}
		if(mergeGreen == 0){
			if(argsP->failDumpB){
				myArgs->failPCC->addThing(anyRes->mainEnt);
				myArgs->failPCC->addThing(anyRes->pairEnt);
			}
			else{
				myArgs->entC->dealloc(anyRes->mainEnt);
				myArgs->entC->dealloc(anyRes->pairEnt);
			}
		}
		myArgs->getPCC->taskCache.dealloc(anyRes);
		anyRes = myArgs->getPCC->getThing();
	}
}

/**Arguments to the merge thread output.*/
typedef struct{
	/**The place to output.*/
	SequenceWriter* curOut;
	/**The place to output, if outputting sam.*/
	CRBSAMFileWriter* samOut;
	/**Arguments to ProSynAr.*/
	ProsynarArgumentParser* argsP;
	/**Output results to write.*/
	ThreadProdComCollector<MergeSequenceData>* resPCC;
	/**Return entry data.*/
	ThreadsafeReusableContainerCache<CRBSAMFileContents>* entC;
} OutputMergeThreadArgs;
/**
 * Output the merged sequences.
 * @param tmpArg The arguments.
 */
void outputMergedResults(void* tmpArg){
	char numBuff[4*sizeof(uintmax_t)+4];
	OutputMergeThreadArgs* myArgs = (OutputMergeThreadArgs*)tmpArg;
	SequenceWriter* curOut = myArgs->curOut;
	CRBSAMFileWriter* samOut = myArgs->samOut;
	MergeSequenceData* anyRes = myArgs->resPCC->getThing();
	while(anyRes){
		if(curOut){
			curOut->nextNameLen = anyRes->seqName.size();
			curOut->nextName = anyRes->seqName.c_str();
			curOut->nextSeqLen = anyRes->seqSeq.size();
			curOut->nextSeq = anyRes->seqSeq.c_str();
			curOut->nextHaveQual = 1;
			curOut->nextQual = &(anyRes->seqQuals[0]);
			curOut->writeNextEntry();
		}
		if(samOut){
			CRBSAMFileContents* curEnt = &(samOut->curEnt);
			curEnt->clear();
			curEnt->lastReadHead = 0;
			curEnt->entryName.insert(curEnt->entryName.end(), anyRes->seqName.begin(), anyRes->seqName.end());
			curEnt->entryFlag = 0;
			curEnt->entryReference.insert(curEnt->entryReference.end(), anyRes->mainEnt->entryReference.begin(), anyRes->mainEnt->entryReference.end());
			curEnt->entryPos = std::min(anyRes->mainEnt->entryPos, anyRes->pairEnt->entryPos);
			curEnt->entryMapq = std::min(anyRes->mainEnt->entryMapq, anyRes->pairEnt->entryMapq);
			sprintf(numBuff, "%ju", (uintmax_t)(curEnt->entrySeq.size()));
			curEnt->entryCigar.insert(curEnt->entryCigar.end(), numBuff, numBuff + strlen(numBuff));
			curEnt->entryCigar.push_back('M');
			curEnt->nextPos = -1;
			curEnt->entryTempLen = 0;
			curEnt->entrySeq.insert(curEnt->entrySeq.end(), anyRes->seqSeq.begin(), anyRes->seqSeq.end());
			curEnt->entryQual.resize(anyRes->seqQuals.size());
			fastaLog10ProbsToPhred(anyRes->seqQuals.size(), &(anyRes->seqQuals[0]), (unsigned char*)(&(curEnt->entryQual[0])));
			samOut->writeNextEntry();
		}
		myArgs->entC->dealloc(anyRes->mainEnt);
		myArgs->entC->dealloc(anyRes->pairEnt);
		myArgs->resPCC->taskCache.dealloc(anyRes);
		anyRes = myArgs->resPCC->getThing();
	}
}

/**Arguments to writing failed sequence.*/
typedef struct{
	/**The place to output.*/
	CRBSAMFileWriter* curOut;
	/**Arguments to ProSynAr.*/
	ProsynarArgumentParser* argsP;
	/**Return entry data.*/
	ThreadsafeReusableContainerCache<CRBSAMFileContents>* entC;
	/**Get tasks to do.*/
	ThreadProdComCollector<CRBSAMFileContents>* taskPCC;
} OutputFailedThreadArgs;
/**
 * Output the failed sequences.
 * @param tmpArg The arguments.
 */
void outputFailedResults(void* tmpArg){
	OutputFailedThreadArgs* myArgs = (OutputFailedThreadArgs*)tmpArg;
	CRBSAMFileWriter* curOut = myArgs->curOut;
	CRBSAMFileContents* anyRes = myArgs->taskPCC->getThing();
	while(anyRes){
		curOut->writeNextEntry(anyRes);
		myArgs->entC->dealloc(anyRes);
		anyRes = myArgs->taskPCC->getThing();
	}
}

#define MAX_QUEUE_SIZE 16

/**
 * Run the damn thing.
 * @param argc The number of arguments.
 * @param argv The arguments.
 * @return Whether there was a problem.
 */
int main(int argc, char** argv){
	int retCode = 0;
	InStream* curInpF = 0;
	TabularReader* curInpT = 0;
	CRBSAMFileReader* curInp = 0;
	OutStream* curOutF = 0;
	SequenceWriter* curOut = 0;
	OutStream* curOutSF = 0;
	TabularWriter* curOutST = 0;
	CRBSAMFileWriter* curOutS = 0;
try{
	//parse arguments, set up
		ProsynarArgumentParser argsP;
			if(argsP.parseArguments(argc-1, argv+1, &std::cout) < 0){
				std::cerr << argsP.argumentError << std::endl;
				retCode = 1;
				goto cleanUp;
			}
			if(argsP.needRun == 0){ goto cleanUp; }
			argsP.performSetup();
	//open the outputs
		if(argsP.seqOutFile){
			openSequenceFileWrite(argsP.seqOutFile, &curOutF, &curOut);
		}
		if(argsP.mergeSamOutFile){
			openCRBSamFileWrite(argsP.mergeSamOutFile, &curOutSF, &curOutST, &curOutS);
		}
	//prepare caches and work queues
		ThreadsafeReusableContainerCache<CRBSAMFileContents> entCache;
		PairedEndCache proCache;
		ThreadProdComCollector<MergeAttemptTask> taskPCC(MAX_QUEUE_SIZE*argsP.numThread);
		ThreadProdComCollector<MergeSequenceData> goodPCC(MAX_QUEUE_SIZE*argsP.numThread);
		ThreadProdComCollector<CRBSAMFileContents> failPCC(MAX_QUEUE_SIZE*argsP.numThread);
	//start up the work threads
		std::vector<MergeThreadArgs> workThrArgs; workThrArgs.resize(argsP.numThread);
		std::vector<void*> liveThread;
		for(intptr_t i = 0; i<argsP.numThread; i++){
			MergeThreadArgs newArg = {(int)i, &argsP, &taskPCC, &entCache, &failPCC, &goodPCC};
			workThrArgs[i] = newArg;
			liveThread.push_back(startThread(attemptMerging, &(workThrArgs[i])));
		}
	//start the good output thread
		OutputMergeThreadArgs goodThrArg = {curOut, curOutS, &argsP, &goodPCC, &entCache};
		void* goodThr = startThread(outputMergedResults, &goodThrArg);
	//start the bad output thread
		CRBSAMFileWriter* failDumpB = argsP.failDumpB;
		OutputFailedThreadArgs failThrArg;
		void* failThr = 0;
		if(failDumpB){
			OutputFailedThreadArgs newFTArg = {failDumpB, &argsP, &entCache, &failPCC};
			failThrArg = newFTArg;
			failThr = startThread(outputFailedResults, &failThrArg);
		}
	//run down the files
		CRBSAMFileContents* curEnt = entCache.alloc();
		uintptr_t numTaskIn = 0;
		bool haveEndHead = false;
		for(uintptr_t si = 0; si < argsP.samNames.size(); si++){
			//open
			openCRBSamFileRead(argsP.samNames[si], &curInpF, &curInpT, &curInp);
			//run down the file looking for unpaired and paired
			while(curInp->readNextEntry(curEnt)){
				//manage the entry
				if(curEnt->lastReadHead){
					if(failDumpB && !haveEndHead){
						failDumpB->writeNextEntry(curEnt);
					}
					if(curOutS && !haveEndHead){
						curOutS->writeNextEntry(curEnt);
					}
				}
				else{
					haveEndHead = true;
					//skip secondary/supplementary stuff
					if(!samEntryIsPrimary(curEnt)){ continue; }
					//if paired, handle
					if(samEntryNeedPair(curEnt)){
						if(proCache.havePair(curEnt)){
							MergeAttemptTask* curPush = taskPCC.taskCache.alloc();
							curPush->mainEnt = curEnt;
							std::pair<uintptr_t,CRBSAMFileContents*> origEnt = proCache.getPair(curEnt);
							curPush->pairEnt = origEnt.second;
							taskPCC.addThing(curPush);
						}
						else{
							proCache.waitForPair(numTaskIn, curEnt);
						}
					}
					else if(failDumpB){
						failPCC.addThing(curEnt);
					}
					else{
						entCache.dealloc(curEnt);
					}
					curEnt = entCache.alloc();
					numTaskIn++;
				}
			}
			//close
			delete(curInp); curInp = 0;
			delete(curInpT); curInpT = 0;
			delete(curInpF); curInpF = 0;
		}
		entCache.dealloc(curEnt); //got one too many
	//drain any outstanding to fail
		std::string badPairName;
		while(proCache.haveOutstanding()){
			std::pair<uintptr_t,CRBSAMFileContents*> origEntP = proCache.getOutstanding();
			CRBSAMFileContents* origEnt = origEntP.second;
			badPairName.clear(); badPairName.insert(badPairName.end(), origEnt->entryName.begin(), origEnt->entryName.end());
			lockMutex(argsP.errLock);
				std::cerr << "Entry " << badPairName << " claims to be paired, but no pair is in file." << std::endl;
			unlockMutex(argsP.errLock);
			if(failDumpB){
				failPCC.addThing(origEnt);
			}
			else{
				entCache.dealloc(origEnt);
			}
		}
	//end the task cache, join the threads
		taskPCC.end();
		for(uintptr_t i = 0; i<liveThread.size(); i++){
			joinThread(liveThread[i]);
		}
	//end the result caches, join the finals
		goodPCC.end();
		failPCC.end();
		joinThread(goodThr);
		if(failThr){ joinThread(failThr); }
}catch(std::exception& err){
	std::cerr << err.what() << std::endl;
	retCode = 1;
}
	cleanUp:
	if(curInp){ delete(curInp); }
	if(curInpT){ delete(curInpT); }
	if(curInpF){ delete(curInpF); }
	if(curOut){ delete(curOut); }
	if(curOutF){ delete(curOutF); }
	if(curOutS){ delete(curOutS); }
	if(curOutST){ delete(curOutST); }
	if(curOutSF){ delete(curOutSF); }
	return retCode;
}

