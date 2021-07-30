#include "prosynar_task.h"

#include <string.h>
#include <algorithm>

#include "whodun_oshook.h"
#include "whodun_parse_seq.h"

ProsynarFilter::ProsynarFilter(){}
ProsynarFilter::~ProsynarFilter(){}

ProsynarMerger::ProsynarMerger(){}
ProsynarMerger::~ProsynarMerger(){}

ProsynarArgumentParser::ProsynarArgumentParser(){
	defOutFN[0] = '-'; defOutFN[1] = 0;
	errLock = makeMutex();
	failDumpS = 0;
	failDumpT = 0;
	failDumpB = 0;
	seqOutFile = 0;
	mergeSamOutFile = 0;
	refFile = 0;
	probFile = 0;
	costFile = 0;
	qualmFile = 0;
	failDumpFile = 0;
	defProbRegMap = 0;
	defAllRegCosts = 0;
	defAllQualMangs = 0;
	numThread = 1;
	useMerger = 0;
	std::map<std::string,ProsynarFilter*(*)()> filtStore;
	getAllProsynarFilters(&filtStore);
	std::map<std::string,ProsynarMerger*(*)()> mergStore;
	getAllProsynarMergers(&mergStore);
	helpDocStore.append("Usage: prosynar [OPTION] [FILE] (-- FILTER)* -- MERGER\n");
	helpDocStore.append("Attempts to merge pairs in SAM files.\n");
	helpDocStore.append("Can filter out questionable pairs, and can select merge algorithms.\n");
	helpDocStore.append("FILTER options are:");
	for(std::map<std::string,ProsynarFilter*(*)()>::iterator tmpIt = filtStore.begin(); tmpIt != filtStore.end(); tmpIt++){
		helpDocStore.push_back(' ');
		helpDocStore.append(tmpIt->first);
	}
	helpDocStore.push_back('\n');
	helpDocStore.append("MERGER options are:");
	for(std::map<std::string,ProsynarMerger*(*)()>::iterator tmpIt = mergStore.begin(); tmpIt != mergStore.end(); tmpIt++){
		helpDocStore.push_back(' ');
		helpDocStore.append(tmpIt->first);
	}
	helpDocStore.push_back('\n');
	helpDocStore.append("The OPTIONS are:\n");
	myMainDoc = helpDocStore.c_str();
	myVersionDoc = "ProSynAr 1.0";
	myCopyrightDoc = "Copyright (C) 2019 UNT HSC Center for Human Identification";
	ArgumentParserStrMeta refMeta("Reference File");
		refMeta.isFile = true;
		refMeta.fileExts.insert(".fa");
		refMeta.fileExts.insert(".fa.gz");
		refMeta.fileExts.insert(".fa.gzip");
		addStringOption("--ref", &refFile, 0, "    Specify the reference sequences in a fasta file.\n    --ref File.fa\n", &refMeta);
	ArgumentParserStrMeta probMeta("Problematic Region File");
		probMeta.isFile = true;
		probMeta.fileExts.insert(".bed");
		addStringOption("--prob", &probFile, 0, "    Specify the problematic regions in a bed file.\n    --prob File.bed\n", &probMeta);
	ArgumentParserStrMeta costMeta("Alignment Parameter File");
		costMeta.isFile = true;
		costMeta.fileExts.insert(".pdc");
		addStringOption("--cost", &costFile, 0, "    Specify the position dependent cost function.\n    --cost File.pdc\n", &costMeta);
	ArgumentParserStrMeta qualmMeta("Quality Mangle File");
		qualmMeta.isFile = true;
		qualmMeta.fileExts.insert(".qualm");
		addStringOption("--qualm", &qualmFile, 0, "    Specify how to modify alignment parameters using their quality.\n    --qualm File.qualm\n", &qualmMeta);
	ArgumentParserStrMeta filDumpMeta("Filtered Dump File");
		filDumpMeta.isFile = true;
		filDumpMeta.fileWrite = true;
		filDumpMeta.fileExts.insert(".sam");
		addStringOption("--faildump", &failDumpFile, 0, "    Specify a file to write reads that were not merged.\n    --faildump File.sam\n", &filDumpMeta);
	ArgumentParserIntMeta threadMeta("Threads");
		addIntegerOption("--thread", &numThread, 0, "    The number of threads to use.\n    --thread 1\n", &threadMeta);
	ArgumentParserStrMeta fastOutMeta("Sequence Output File");
		fastOutMeta.isFile = true;
		fastOutMeta.fileWrite = true;
		fastOutMeta.fileExts.insert(".fq");
		fastOutMeta.fileExts.insert(".fq.gz");
		fastOutMeta.fileExts.insert(".fq.gzip");
		fastOutMeta.fileExts.insert(".fastq");
		fastOutMeta.fileExts.insert(".fastq.gz");
		fastOutMeta.fileExts.insert(".fastq.gzip");
		addStringOption("--out", &seqOutFile, 0, "    Specify a raw sequence output file.\n    --out File.fq\n", &fastOutMeta);
	ArgumentParserStrMeta samOutMeta("Alignment Hint Output File");
		samOutMeta.isFile = true;
		samOutMeta.fileWrite = true;
		samOutMeta.fileExts.insert(".sam");
		addStringOption("--samze", &mergeSamOutFile, 0, "    Specify a location to write merged alignments.\n    These are zero effort alignments: they are simply placement in the genome.\n    These will need realignment.\n    --out File.sam\n", &samOutMeta);
}

ProsynarArgumentParser::~ProsynarArgumentParser(){
	if(errLock){ killMutex(errLock); }
	for(uintptr_t i = 0; i<useFilters.size(); i++){ delete(useFilters[i]); }
	if(useMerger){ delete(useMerger); }
	if(failDumpB){ delete(failDumpB); }
	if(failDumpT){ delete(failDumpT); }
	if(failDumpS){ delete(failDumpS); }
}

int ProsynarArgumentParser::handleUnknownArgument(int argc, char** argv, std::ostream* helpOut){
	if(strcmp(argv[0],"--")==0){
		std::map<std::string,ProsynarFilter*(*)()> filtStore;
		getAllProsynarFilters(&filtStore);
		std::map<std::string,ProsynarMerger*(*)()> mergStore;
		getAllProsynarMergers(&mergStore);
		int workAC = argc - 1;
		char** workAV = argv + 1;
		while(workAC){
			if(strcmp(*workAV,"--")==0){
				workAC--; workAV++;
				continue;
			}
			std::map<std::string,ProsynarFilter*(*)()>::iterator filtIt = filtStore.find(workAV[0]);
			std::map<std::string,ProsynarMerger*(*)()>::iterator mergIt = mergStore.find(workAV[0]);
			workAC--; workAV++;
			if(filtIt != filtStore.end()){
				ProsynarFilter* newFilt = (filtIt->second)();
				useFilters.push_back(newFilt);
				int numAte = newFilt->parseArguments(workAC, workAV, helpOut);
				if(numAte < 0){ argumentError = newFilt->argumentError; }
				else{ workAC -= numAte; workAV += numAte; }
				if((numAte < 0) || (newFilt->needRun == 0)){ needRun = 0; return (numAte < 0) ? -1 : (workAV - argv); }
			}
			else if(mergIt != mergStore.end()){
				if(useMerger){
					argumentError = "Multiple mergers requested.";
					return -1;
				}
				useMerger = (mergIt->second)();
				int numAte = useMerger->parseArguments(workAC, workAV, helpOut);
				if(numAte < 0){ argumentError = useMerger->argumentError; }
				else{ workAC -= numAte; workAV += numAte; }
				if((numAte < 0) || (useMerger->needRun == 0)){ needRun = 0; return (numAte < 0) ? -1 : (workAV - argv); }
			}
			else{
				argumentError = "Unknown filter/merger: ";
				argumentError.append(workAV[0]);
				return -1;
			}
		}
		return workAV - argv;
	}
	else if(strcmp(argv[0], "--helplistm")==0){
		needRun = 0;
		std::map<std::string,ProsynarMerger*(*)()> mergStore;
		getAllProsynarMergers(&mergStore);
		for(std::map<std::string,ProsynarMerger*(*)()>::iterator mrgIt = mergStore.begin(); mrgIt != mergStore.end(); mrgIt++){
			(*helpOut) << mrgIt->first << std::endl;
		}
		return 1;
	}
	else if(strcmp(argv[0], "--helplistf")==0){
		needRun = 0;
		std::map<std::string,ProsynarFilter*(*)()> filtStore;
		getAllProsynarFilters(&filtStore);
		for(std::map<std::string,ProsynarFilter*(*)()>::iterator mrgIt = filtStore.begin(); mrgIt != filtStore.end(); mrgIt++){
			(*helpOut) << mrgIt->first << std::endl;
		}
		return 1;
	}
	else{
		samNames.push_back(argv[0]);
		return 1;
	}
}
void ProsynarArgumentParser::printExtraGUIInformation(std::ostream* toPrint){
	(*toPrint) << "STRINGVEC\t<|>\tNAME\tInput Files\tFILE\tREAD\t3\t.sam\t.sam.gz\t.sam.gzip\tDOCUMENT\t53414D2066696C657320746F207265616C69676E2E" << std::endl;
}
int ProsynarArgumentParser::posteriorCheck(){
	if(numThread <= 0){
		argumentError = "thread must be positive.";
		return 1;
	}
	if(useMerger == 0){
		argumentError = "No merge operation specified.";
		return 1;
	}
	if(useMerger->posteriorCheck()){
		argumentError = useMerger->argumentError;
		return 1;
	}
	for(uintptr_t i = 0; i<useFilters.size(); i++){
		if(useFilters[i]->posteriorCheck()){
			argumentError = useFilters[i]->argumentError;
			return 1;
		}
	}
	if(samNames.size() == 0){
		samNames.push_back("-");
	}
	if(seqOutFile && (strlen(seqOutFile)==0)){
		seqOutFile = 0;
	}
	if(failDumpFile && (strlen(failDumpFile)==0)){
		failDumpFile = 0;
	}
	if(mergeSamOutFile && (strlen(mergeSamOutFile)==0)){
		mergeSamOutFile = 0;
	}
	if(!seqOutFile && !mergeSamOutFile){
		seqOutFile = defOutFN;
	}
	return 0;
}

void ProsynarArgumentParser::performSetup(){
	std::string fileConts;
	//load the reference
	if(refFile){
		InStream* redFile = 0;
		SequenceReader* refRead = 0;
		try{
			//reference
				openSequenceFileRead(refFile, &redFile, &refRead);
				while(refRead->readNextEntry()){
					std::string crefNam(refRead->lastReadName, refRead->lastReadName + refRead->lastReadShortNameLen);
					std::string crefSeq(refRead->lastReadSeq, refRead->lastReadSeq + refRead->lastReadSeqLen);
					allRefs[crefNam] = crefSeq;
				}
				delete(refRead); refRead = 0;
				delete(redFile); redFile = 0;
		}catch(...){
			if(refRead){ delete(refRead); }
			if(redFile){ delete(redFile); }
			throw;
		}
	}
	//load in the costs
	if(costFile){
		InStream* redFile = new FileInStream(costFile);
		fileConts.clear(); readStream(redFile, &fileConts);
		delete(redFile); redFile = 0;
		parseMultiregionPositionDependentCost(fileConts.c_str(), fileConts.c_str() + fileConts.size(), &allRegCosts);
		defAllRegCosts = &allRegCosts;
	}
	//load in the quality mangles
	if(qualmFile){
		InStream* redFile = new FileInStream(qualmFile);
		fileConts.clear(); readStream(redFile, &fileConts);
		delete(redFile); redFile = 0;
		parseMultiregionPositionQualityMangle(fileConts.c_str(), fileConts.c_str() + fileConts.size(), &allQualMangs);
		defAllQualMangs = &allQualMangs;
	}
	//problem regions
	if(probFile){
		InStream* redFile = 0;
		TabularReader* tabRead = 0;
		try{
			redFile = new FileInStream(probFile);
			tabRead = new TSVTabularReader(0, redFile);
			BedFileReader proBed(tabRead);
			delete(tabRead); tabRead = 0;
			delete(redFile); redFile = 0;
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
		}catch(...){
			if(redFile){ delete(redFile); }
			if(tabRead){ delete(tabRead); }
			throw;
		}
		defProbRegMap = &probRegMap;
	}
	//let the filters and merger prepare
	for(uintptr_t i = 0; i<useFilters.size(); i++){
		useFilters[i]->initialize(this);
	}
	useMerger->initialize(this);
	//open up the fail dump, if any
	if(failDumpFile){
		openCRBSamFileWrite(failDumpFile, &failDumpS, &failDumpT, &failDumpB);
	}
}

//*****************************************************************************
//Add new filters/merge algorithms here.

#include "prosynar_task_pear.h"
#include "prosynar_task_flash.h"
#include "prosynar_task_naimerge.h"
#include "prosynar_task_refover.h"
#include "prosynar_task_refprob.h"

void getAllProsynarFilters(std::map<std::string,ProsynarFilter*(*)()>* toFill){
	(*toFill)["Frover"] = factoryReferenceOverlapFilter;
	(*toFill)["Fprover"] = factoryProbabilisticReferenceOverlapFilter;
	(*toFill)["Fpprobreg"] = factoryProblematicRegionFilter;
}

void getAllProsynarMergers(std::map<std::string,ProsynarMerger*(*)()>* toFill){
	(*toFill)["Mpear"] = factoryPearMerger;
	(*toFill)["Mflash"] = factoryFlashMerger;
	(*toFill)["Malign"] = factorySimpleAlignMerger;
}
