

class NegCtrlClass{
public:
	FeatureStruct *negctrls;
int NofNegCtrls;
void InitialiseData(void);
void FillNegativeCtrls();
void WriteBinCoverage(string,int);
void ReadBinCoverage(string,int); // If the SAMFILE already processed, read the number of reads per bin
void AssociateProbeswithNegativeControls(ProbeSet&);
void BinPeaks(int,int,int);

private:
	void BinInteractor(int,int,int);

};


void NegCtrlClass::InitialiseData(){
	int i=0;
	NofNegCtrls=366;
	negctrls = new FeatureStruct [400];
}

void NegCtrlClass::FillNegativeCtrls(){

#ifdef UNIX
string filename1,filename2;
string dirname="/bubo/proj/b2011029/bin/3CAnalysis/";
filename1.append(dirname);
filename2.append(dirname);
filename1.append("100exons_min100kfromTSS_GATC_150perend.bed");
ifstream infile1(filename1.c_str());
filename2.append("100intergenic_min100kfromTSS_GATC_150perend.bed");
ifstream infile2(filename2.c_str());
#endif

#ifdef WINDOWS
	ifstream infile1("100exons_min100kfromTSS_GATC_150perend.bed");
	ifstream infile2("100intergenic_min100kfromTSS_GATC_150perend.bed");
#endif
int i=0;
string chr;
do{
	infile1 >>	negctrls[i].chr >> negctrls[i].start >> negctrls[i].end;
	negctrls[i].midpoint=(abs((negctrls[i].end-negctrls[i].start))/2)+negctrls[i].start;
	negctrls[i].type="G";
	++i;
}while(!infile1.eof());

infile1.close();
i--;

do{
	infile2 >> negctrls[i].chr >> negctrls[i].start >> negctrls[i].end;
	negctrls[i].midpoint=(abs((negctrls[i].end-negctrls[i].start))/2)+negctrls[i].start;
	negctrls[i].type="I";
	++i;

}while(!infile2.eof());
infile2.close();
NofNegCtrls=i-1;
/*
	ifstream infile1("Probes_notAssociatedwithAnyFeatures.txt");

int i=0;
string chr;
do{
	infile1 >>	negctrls[i].chr >> negctrls[i].start >> negctrls[i].end;
	negctrls[i].midpoint=(abs((negctrls[i].end-negctrls[i].start))/2)+negctrls[i].start;
	negctrls[i].type="G";
	++i;
}while(!infile1.eof());

infile1.close();
i--;

NofNegCtrls=i-1;
*/

	for(i=0;i<NofNegCtrls;++i){
		negctrls[i].AllPeaks_PeakBins.resize(NumberofPeakFiles);
		negctrls[i].AllExperiments_IntBins.resize(NumberofExperiments);
		for(int e=0; e<NumberofExperiments;++e){
			negctrls[i].AllExperiments_IntBins[e].interactorbins.push_back(vector <int> ());
			negctrls[i].AllExperiments_IntBins[e].Interactor.push_back(vector <double>());
			negctrls[i].AllExperiments_IntBins[e].EnrichedBins.push_back(vector <int> ());
		}
		for(int j=0;j<NumberofPeakFiles; ++j){
			negctrls[i].AllPeaks_PeakBins[j].PeakBins.push_back(vector <int> ());
			negctrls[i].AllPeaks_PeakBins[j].PeakBins[0].resize(NumberofBins);
		}
	}
	for(i=0;i<NofNegCtrls;++i){
		for(int e=0; e<NumberofExperiments;++e){
			negctrls[i].AllExperiments_IntBins[e].interactorbins[0].resize(NumberofBins);
			negctrls[i].AllExperiments_IntBins[e].Interactor[0].resize(NumberofBins);
		}
	}
}

void NegCtrlClass::BinInteractor(int negctrlid, int interactorcoord,int whichExperiment){
	
	int diff,bin;
	diff=(interactorcoord-negctrls[negctrlid].midpoint);
	if(abs(diff)<MaxJunctionDistance){		
		bin=((diff/BinSize)+NofInteractorBins);
		negctrls[negctrlid].AllExperiments_IntBins[whichExperiment].interactorbins[0][bin-1]++;
	}
}
void NegCtrlClass::BinPeaks(int negctrlid, int interactorcoord,int whichAb){
	int diff,bin;
	diff=(interactorcoord-negctrls[negctrlid].midpoint);
	if(abs(diff)<MaxJunctionDistance){		
		bin=((diff/BinSize)+NofInteractorBins);
		negctrls[negctrlid].AllPeaks_PeakBins[whichAb].PeakBins[0][bin]++;
	}
}


void NegCtrlClass::AssociateProbeswithNegativeControls(ProbeSet& mm9prs){


int i=0,j=0,k=0,x=0;
int diff;
long int startsearch,endsearch;


for(i=0;i<NofNegCtrls;++i){
	startsearch=0;endsearch=-1;
	for(j=0;j<mm9prs.ChrNames.size();++j){
		if(negctrls[i].chr==mm9prs.ChrNames[j]){
			startsearch=mm9prs.ChrRowStartIndexes[j];
			endsearch=mm9prs.ChrRowEndIndexes[j];
			break;
		}
	}
	x=0;
	for(j=startsearch;j<=endsearch;++j){
		diff=mm9prs.MM9Probes[j].center-negctrls[i].midpoint;
		if(abs(diff)<AssociateInteractions){
			if(mm9prs.MM9Probes[j].annotated==0){
				mm9prs.MM9Probes[j].annotated=3; // Associate to a negative control
				negctrls[i].ProbeIDs.push_back(j);
			}
			else
				mm9prs.MM9Probes[j].conflicting_annotations=1;
		}
		if(diff>BinSize)
			break;
	}
	if(i%100==0)
		cout << i << "    Probes Associated with Negative Controls" << endl;
}
}

void NegCtrlClass::ReadBinCoverage(string InputFileNameBase,int whichexperiment){ // If the SAMFILE already processed, read the number of reads per bin
	int i,j;
	string infilename;
	
	infilename.append(InputFileNameBase);
	infilename.append("CoveragePerBin_negctrls.txt");
	ifstream infile(infilename.c_str());
	for(i=0;i<NofNegCtrls;++i){
		for(j=0;j<NumberofBins;++j)
			infile >> negctrls[i].AllExperiments_IntBins[whichexperiment].interactorbins[0][j];
	}
	infile.close();
}
void NegCtrlClass::WriteBinCoverage(string OutputFileNameBase,int whichexperiment){
	int i,j;
	string outfilename;
	
	outfilename.append(OutputFileNameBase);
	outfilename.append("CoveragePerBin_negctrls.txt");
	ofstream outfile(outfilename.c_str());
	
	for(i=0;i<NofNegCtrls;++i){
		for(j=0;j<NumberofBins;++j)
			outfile << negctrls[i].AllExperiments_IntBins[whichexperiment].interactorbins[0][j] << '\t';
		outfile << endl;
	}
	outfile.close();
}
