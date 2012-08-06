/* 
    DetectInteractions - Detect Genome Interactions
    Copyright (C) 2012  Pelin Akan (pelin.akan@scilifelab.se)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.                  

    For a copy of the GNU Affero General Public License
    see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cstring>
#include <string>
#include <time.h>
#include <vector>
#include <omp.h>
#include <sstream>
using namespace std;
#include "ap.h"

//#define UNIX
#define WINDOWS

#define SAM_FILE_ALREADY_PROCESSED // Only read the interactor matrix for features
//#define PROCESS_SAMFILE //Process the SAMFILE
//#define PROCESS_SAM_FILES // Process all samfiles to generate interactor matrices for different bin sizes 

// disable some irrelevant warnings
#if (AE_COMPILER==AE_MSVC)
#pragma warning(disable:4100)
#pragma warning(disable:4127)
#pragma warning(disable:4702)
#pragma warning(disable:4996)
#endif
#include "alglibmisc.h"
#include "alglibinternal.h"
#include "linalg.h"
#include "statistics.h"
#include "dataanalysis.h"
#include "specialfunctions.h"
#include "solvers.h"
#include "optimization.h"
#include "diffequations.h"
#include "fasttransforms.h"
#include "integration.h"
#include "interpolation.h"
using namespace alglib_impl;

const int NUM_OF_THREADS=8;
const int BUFFERSIZE=50000;
const int NumberofPeakFiles=7;
const double SignificanceThreshold=0.001;
const int MinNumberofPairs=2;
const int MinNumberofPeaksPerBin=0;
const int NumberofGenes=24900;
const int NofNegCtrls=400;
const int NofIsoforms=20;
const int ClusterPromoters=5000; //Cluster Promoters of Isoforms that are 5 kb away from each other
const int AssociateInteractions=5000; //decide if an interaction comes from a feature or not
const int ExcludeOutliers=10;
const int NumberofExperiments = 1;
int NofInteractorBins; // Depends on the maximum allowed junction distance
int MaxJunctionDistance=100000;
int BinSize=5000;
int NumberofBins;
double ExpressionThr=2.0;
int CellType=0; // 0:mES, 1:XEN, 2:TS
int padding=700; //For Sequence Capture Probes


#include "DataStructs.h"
#include "AssociateProbeswithFeatures.h"
#include "PromoterClass.h"
#include "NegativeControls.h"
#include "GetOptions.h"
#include "AnnotateProms_with_PeakFiles.h"
#include "GenerateFileNames.h"
#include "linear.h"
#include "BackgroundInteractionFrequency.h"
#include "DetectEnrichedBins.h"
#include "VisualiseInteractions.h"

void ProcessPeaks(PromoterClass& Promoters, NegCtrlClass& NegativeControls,PeakClass& PC,DetectEnrichedBins& EnrichedBins,string &PeakFN, string abname, bool ifBED,string INTERACTIONFILENAMEBASE,int CellType,int abindex,int ExperimentIndex){
	ifstream PeakFile;
	PeakFile.open(PeakFN.c_str());
	PC.ReadPeakFile(PeakFile,ifBED);
	PeakFile.close();
	
	string FileName1, FileName2,FileName3,FileName4,FileName5;
	
	PC.AnnotatewithNegCtrls(NegativeControls,abindex);
	PC.AnnotatewithPromoters(Promoters,abindex);
	EnrichedBins.AssociatePeaksWithIntBins(Promoters,INTERACTIONFILENAMEBASE,abname,CellType,abindex,ExperimentIndex);
	EnrichedBins.AssociatePeaksWithIntBins_NegCtrls(NegativeControls,INTERACTIONFILENAMEBASE,abname,abindex,ExperimentIndex);

}

int main (int argc,char* argv[]){
string INTERACTIONFILENAMEBASE, BaseFileName;
#ifdef UNIX
 //FOR UNIX

	if (argc < 5) {
		print_usage();
		return -1;
	}
	INTERACTIONFILENAMEBASE.append(argv[1]);
	CellType=atoi(argv[2]);
	BinSize=atoi(argv[3]);
	MaxJunctionDistance=(atoi(argv[4]));
	NofInteractorBins=(MaxJunctionDistance)/BinSize;
	NumberofBins=NofInteractorBins*2;
	cout << "INTERACTORFILENAMEBASE" << INTERACTIONFILENAMEBASE << endl;
	cout << "Cell Type        " <<  CellType << endl;
	cout << "Bin Size         " << BinSize     << endl;
	cout << "Max Junction Distance " << MaxJunctionDistance << endl;
	cout << "Number of Bins on Each Side    " << NofInteractorBins << endl;


#endif
#ifdef WINDOWS
//	INTERACTIONFILENAMEBASE="3C-mES_BR1_5kbBins_MaxJD100kb_";
	BinSize=5000;
	MaxJunctionDistance=500000;
	NofInteractorBins=(MaxJunctionDistance)/BinSize;
	NumberofBins=NofInteractorBins*2;
#endif
//   --        INITIALISE CLASSES   --
	PromoterClass Promoters;
	NegCtrlClass NegativeControls;
	ProbeSet mm9probes;
	DetermineBackgroundLevels BackgroundLevels;
	DetectEnrichedBins EnrichedBins;
	vector < PeakClass > AllPeaks;
	VisualiseInteractions UCSCTracks;
//-------------------//------------------------------
string ext1,ext2,FileName1,FileName2,FileName3,FileName4;

//Read Promoters and Negative Controls, Annotate
	Promoters.InitialiseData();
	Promoters.ReadPromoterAnnotation();
	NegativeControls.InitialiseData();
	NegativeControls.FillNegativeCtrls();
	cout << "NegativeControls Read" << endl;

// ASSOCIATE PROBES WITH FEATURES
	mm9probes.ReadProbeCoordinates();
	mm9probes.InitialiseData();
	Promoters.AssociateProbeswithPromoters(mm9probes);
	NegativeControls.AssociateProbeswithNegativeControls(mm9probes);
//	SAMFILE.open(SAMFILENAME.c_str());
//	samfile.ProcessTheSAMFile_NoCTX_Probes(mm9probes);
//	SAMFILE.close();
//	mm9probes.PrinttheBins(BaseFileName);


#ifdef SAM_FILE_ALREADY_PROCESSED
	ifstream ExperimentsFile("Experiments.txt");
	for(int experimentindex=0; experimentindex<NumberofExperiments;++experimentindex){
		ExperimentsFile >> INTERACTIONFILENAMEBASE;
		Promoters.ReadBinCoverage(INTERACTIONFILENAMEBASE,experimentindex);
		NegativeControls.ReadBinCoverage(INTERACTIONFILENAMEBASE,experimentindex);
		//CALCULATE BACKGROUND LEVELS	
		cout << INTERACTIONFILENAMEBASE << "       Background Levels Calculated" << endl;		
		BackgroundLevels.InitialiseVars();
		BackgroundLevels.CalculateMeanandStd(Promoters, NegativeControls,INTERACTIONFILENAMEBASE,experimentindex);
		//DETECT INTERACTIONS ABOVE BACKGROUND LEVELS FOR EACH BIN
		EnrichedBins.InitialiseData();
		EnrichedBins.DetectEnrichedInteractionBins(Promoters,BackgroundLevels.mean,BackgroundLevels.stdev,INTERACTIONFILENAMEBASE,experimentindex);
		EnrichedBins.DetectEnrichedInteractionBins_NegCtrls(NegativeControls,BackgroundLevels.mean,BackgroundLevels.stdev,INTERACTIONFILENAMEBASE,experimentindex);
		cout << INTERACTIONFILENAMEBASE <<  "      Enriched Bins Detected" << endl;
		//READ, PROCESS THE PEAKS
		string PeakFileName;
		int abindex=0;
		ifstream AbNameFile;
		AbNameFile.open("Abnames.txt");
		do{
			AbNameFile >> PeakFileName >> ext1;
			if(PeakFileName.compare("END")==0)
				break;
#ifdef UNIX
			PeakFileName.append("/bubo/proj/b2011029/bin/3CAnalysis/");
#endif	
			AbNames.push_back(ext1);
			PeakClass Peaks;
			ProcessPeaks(Promoters,NegativeControls,Peaks,EnrichedBins,PeakFileName,ext1,0,INTERACTIONFILENAMEBASE,CellType,abindex,experimentindex);
			AllPeaks.push_back(Peaks);
			Peaks.~PeakClass();
			AllPeaks[abindex].abnames.push_back(ext1);
			++abindex;
			PeakFileName.clear();
			ext1.clear();
		}while(PeakFileName!="END");
		cout << INTERACTIONFILENAMEBASE << "          All Peak Files Read" << endl;
		
		EnrichedBins.PrintMetaAssociationwithPeaks(Promoters,NegativeControls,INTERACTIONFILENAMEBASE,experimentindex);

		AbNameFile.close();
		BackgroundLevels.~DetermineBackgroundLevels();
		EnrichedBins.~DetectEnrichedBins();
		INTERACTIONFILENAMEBASE.clear(); 
	}
#ifdef UNIX		
		ofstream GFFFILE("/bubo/proj/b2011029/bin/3CAnalysis/GFF_File.txt");
#endif
		ofstream GFFFile("GFF_file.txt");
		UCSCTracks.WriteGFFFiles(Promoters,GFFFile);
		ofstream BEDFile("BED_file.txt");
#ifdef UNIX		
		ofstream BEDFILE("/bubo/proj/b2011029/bin/3CAnalysis/BED_file.txt");
#endif
		UCSCTracks.WriteBEDFiles(Promoters,BEDFile);
#endif
}