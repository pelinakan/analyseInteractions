struct temppars{
		string chr;
		int start;
		int end;
		double expression[3];
		string strand;
		string name;
		string tr_id;
};
class PromoterClass{ //Probe Clusters Associated with a Promoter
	friend class NegCtrlClass;
	friend class ProbeSet;
	friend class RESitesClass;
	vector <string> ChrNames_refseq;
	vector < vector <int> > refseq_indexes; //Based on refseq gene chromosome
public:
FeatureStruct *refseq;	
temppars *tp;
void InitialiseData(void);
void ReadPromoterAnnotation();
void WriteBinCoverage(string,int);
void ReadBinCoverage(string, int, RESitesClass&, MappabilityClass&); // If the SAMFILE already processed, read the number of reads per bin
void BinPromoterInteractor(int,int,int,int);
void BinPeaks(int,int,int,int);
void PrintEnrichments_PerTSS(string,int);
void AssociateProbeswithPromoters(ProbeSet&);
bool AnnotateWithPromoters(string,int);
int NofPromoters;
private:
	void GetTrFeats(stringstream&,temppars&);
	void ClusterIsoformPromoters(vector<int>,vector<string>,int);
	void GetPairInformation(stringstream&, stringstream&, PairStruct&);
};

void PromoterClass::InitialiseData(void){

	tp=new temppars [2];
	refseq = new FeatureStruct [NumberofGenes];
	for(int i=0;i<NumberofGenes;i++){
		refseq[i].expression=(double *) calloc(3,sizeof(double));
	}	
}
void PromoterClass::GetTrFeats(stringstream &trx, temppars &tpars){

	string field;
	getline(trx,tpars.name,'\t');

	getline(trx,tpars.tr_id,'\t'); // tr id
	getline(trx,tpars.chr,'\t');
	getline(trx,tpars.strand,'\t');
	getline(trx,field,'\t');
	if(tpars.strand=="+"){
		tpars.start=atoi(field.c_str());
		getline(trx,field,'\t');
		tpars.end=atoi(field.c_str());
	}
	else{
		tpars.end=atoi(field.c_str());
		getline(trx,field,'\t');
		tpars.start=atoi(field.c_str());
		}
	getline(trx,field,'\t');
	getline(trx,field,'\t');
	getline(trx,field,'\t');

	getline(trx,field,'\t');
	tpars.expression[0]=atof(field.c_str()); //mES
	getline(trx,field,'\t');
	tpars.expression[1]=atof(field.c_str()); //XEN
	getline(trx,field,'\t');
	tpars.expression[2]=atof(field.c_str()); //TS
	

}
void PromoterClass::ClusterIsoformPromoters(vector <int> isoformprs, vector<string> tr_ids, int gene_idx){
unsigned int j,k,l,cluster_idx=0;
vector<int> clustercoords;
vector<string> clustered_tr_ids;
vector<bool> clustered;

for(j=0;j<isoformprs.size();++j)
	clustered.push_back(0);

for(j=0;j<isoformprs.size();++j){
	l=0;
	if(!clustered[j]){
		clustercoords.push_back(isoformprs[j]);
		clustered_tr_ids.push_back(tr_ids[j]);
		if(j+1<isoformprs.size()){
			for(k=j+1;k<isoformprs.size();++k){
				if(!clustered[j]){
					if(abs(isoformprs[j]-isoformprs[k])<ClusterPromoters){
						++l;
						clustered[k]=1;
						clustercoords.push_back(isoformprs[k]);
						clustered_tr_ids.push_back(tr_ids[k]);
					}
				}
			}
			refseq[gene_idx].isoformpromotercoords.push_back(clustercoords[0]+((clustercoords[l]-clustercoords[0])/2));
			refseq[gene_idx].TranscriptName.push_back(clustered_tr_ids[0]);
			++cluster_idx;
			clustercoords.clear();
			clustered_tr_ids.clear();
		}
		else{
			refseq[gene_idx].isoformpromotercoords.push_back(isoformprs[j]);
			refseq[gene_idx].TranscriptName.push_back(tr_ids[j]);
		}
	}
}
}

void PromoterClass::ReadPromoterAnnotation()
{
	string temp,tr1,tr2;
	int geneindex=0;

#ifdef UNIX
string RefSeqfilename;
 RefSeqfilename.append("/bubo/proj/b2011029/bin/Detect_Interactions/supportingFiles/");
RefSeqfilename.append("mm9_refseq_sortedbyname_withexpression.txt");
ifstream RefSeq_file(RefSeqfilename.c_str());
#endif
#ifdef WINDOWS
	ifstream RefSeq_file("mm9_refseq_NMonly_sortedbyname_withexpression.txt");
#endif
	vector <int> isoformprs; // to keep isoform promoters
	vector <string> tr_ids;
	getline(RefSeq_file,temp); //get the header row
	
	getline(RefSeq_file,tr1);
	stringstream trx1 ( tr1 ); 
	GetTrFeats(trx1,tp[0]);

	isoformprs.push_back(tp[0].start);
	tr_ids.push_back(tp[0].tr_id);
	do{
		getline(RefSeq_file,tr2);
		stringstream trx2 ( tr2);
		GetTrFeats(trx2,tp[1]);
		while(tp[0].name==tp[1].name){
			isoformprs.push_back(tp[1].start);
			tr_ids.push_back(tp[1].tr_id);
			getline(RefSeq_file,tr2);
			stringstream trx2 ( tr2);
			GetTrFeats(trx2,tp[1]);
		};
		if(isoformprs.size()>1)
			ClusterIsoformPromoters(isoformprs,tr_ids,geneindex); //Promoters that are within 5 kb of each other are clustered as one
		else{
			refseq[geneindex].isoformpromotercoords.push_back(isoformprs[0]);
			refseq[geneindex].TranscriptName.push_back(tp[0].tr_id);
		}
		refseq[geneindex].RefSeqName.append(tp[0].name);
		refseq[geneindex].chr.append(tp[0].chr);
		refseq[geneindex].gene_end=tp[0].end;
		refseq[geneindex].strand=tp[0].strand;
		refseq[geneindex].expression[0]=tp[0].expression[0];
		refseq[geneindex].expression[1]=tp[0].expression[1];
		refseq[geneindex].expression[2]=tp[0].expression[2];
		++geneindex;
		isoformprs.clear();
		isoformprs.push_back(tp[1].start);
		tr_ids.clear();
		tr_ids.push_back(tp[1].tr_id);
//SWAP TP[0] and TP[1]
		tp[0].name.clear();
		tp[0].name.append(tp[1].name);
		tp[0].tr_id.clear();
		tp[0].tr_id.append(tp[1].tr_id);
		tp[0].chr.clear();
		tp[0].chr.append(tp[1].chr);
		tp[0].end=tp[1].end;
		tp[0].start=tp[1].start;
		tp[0].strand=tp[1].strand;
		tp[0].expression[0]=tp[1].expression[0];
		tp[0].expression[1]=tp[1].expression[1];
		tp[0].expression[2]=tp[1].expression[2];

		if(geneindex%10000==0)
		  //break;
		cout << geneindex << "    Promoters Annotated" << endl;

	}while(tp[1].name!="END");

	NofPromoters=geneindex;

	// Index Promoters for faster access
	unsigned int i=0;
	int found=0;
	for(geneindex=0;geneindex<NofPromoters;++geneindex){
		found=-1;
		for(i=0;i<ChrNames_refseq.size();++i){
			if(refseq[geneindex].chr.compare(ChrNames_refseq[i])==0)
				found=i;
		}
		if(found==-1){
			ChrNames_refseq.push_back(refseq[geneindex].chr);
			refseq_indexes.push_back(vector<int>());
			refseq_indexes[ChrNames_refseq.size()-1].push_back(geneindex);
		}
		else
			refseq_indexes[found].push_back(geneindex);

	}
	
	cout << "RefSeq Chromosome Indexes Generated for   " << NofPromoters <<  endl;

	//Initialise Promoter Feature vectors
	for(i=0;i<NofPromoters;++i){
		refseq[i].AllPeaks_PeakBins.resize(NumberofPeakFiles);
		refseq[i].AllExperiments_IntBins.resize(NumberofExperiments);
		for(int e=0;e<NumberofExperiments;++e){
			for(int j=0;j<refseq[i].isoformpromotercoords.size();++j){
				refseq[i].AllExperiments_IntBins[e].interactorbins.push_back(vector <double> ());
				refseq[i].AllExperiments_IntBins[e].Interactor.push_back(vector <double>());
				refseq[i].AllExperiments_IntBins[e].ClusteredEnrichedBins.push_back(vector <double> ());
			}
		}
		for(int k=0;k<NumberofPeakFiles;++k)
			for(int j=0;j<refseq[i].isoformpromotercoords.size();++j)
				refseq[i].AllPeaks_PeakBins[k].PeakBins.push_back(vector <int> ());
	}
	for(i=0;i<NofPromoters;++i){
		for(int e=0;e<NumberofExperiments;++e){
			for(int j=0;j<refseq[i].isoformpromotercoords.size();++j){
				refseq[i].AllExperiments_IntBins[e].interactorbins[j].resize(NumberofBins);
				refseq[i].AllExperiments_IntBins[e].Interactor[j].resize(NumberofBins);
			}
		}
		for(int k=0;k<NumberofPeakFiles;++k)
			for(int j=0;j<refseq[i].isoformpromotercoords.size();++j)
				refseq[i].AllPeaks_PeakBins[k].PeakBins[j].resize(NumberofBins);
	}
}



void PromoterClass::BinPromoterInteractor(int promid, int interactorcoord, int whichisoform, int whichexperiment){

	int diff,bin;
	diff=(interactorcoord-refseq[promid].isoformpromotercoords[whichisoform]);
	if(abs(diff)<MaxJunctionDistance){		
		if(refseq[promid].strand=="+"){
			bin=((diff/BinSize)+NofInteractorBins);
			refseq[promid].AllExperiments_IntBins[whichexperiment].interactorbins[whichisoform][bin]++;
		}
		else{
			bin=((NofInteractorBins-diff/BinSize));		
			refseq[promid].AllExperiments_IntBins[whichexperiment].interactorbins[whichisoform][bin]++;
		}
	}

}
void PromoterClass::BinPeaks(int promid, int interactorcoord, int whichisoform,int whichAb){

	int diff,bin;
	diff=(interactorcoord-refseq[promid].isoformpromotercoords[whichisoform]);
	if(refseq[promid].strand=="+"){
		bin=((diff/BinSize)+NofInteractorBins);
		refseq[promid].AllPeaks_PeakBins[whichAb].PeakBins[whichisoform][bin]++;
	}
	else{
		bin=((NofInteractorBins-diff/BinSize));		
		refseq[promid].AllPeaks_PeakBins[whichAb].PeakBins[whichisoform][bin]++;
	}
}

void PromoterClass::PrintEnrichments_PerTSS(string BaseFileName,int CellType){

	int i,j,nofabs;
	string fn1, fn2;
	fn1.append(BaseFileName);
	fn1.append("_EnrichmentperTSS_expressed.txt");
	fn2.append(BaseFileName);
	fn2.append("_EnrichmentperTSS_notexpressed.txt");
	
	ofstream outf1(fn1.c_str());
	ofstream outf2(fn2.c_str());

	nofabs=refseq[0].EnrPerTSS.peakbins.size(); //Nof Antibodies to probe possible enhances

	for(i=0;i<NofPromoters;++i){
//		outf1 << refseq[i].RefSeqName << '\t';
		if(refseq[i].expression[CellType]>ExpressionThr){
			for(j=0;j<nofabs;++j){
				outf1 << refseq[i].EnrPerTSS.peakbins[j] << '\t';
				outf1 << refseq[i].EnrPerTSS.nofpeakbins[j] << '\t';
				outf1 << refseq[i].EnrPerTSS.nonpeakbins[j] << '\t';
				outf1 << refseq[i].EnrPerTSS.nofnonpeakbins[j] << '\t';
			}
			outf1 << endl;
		}
		else{
			if(refseq[i].expression[CellType]!=1.00E-20){
				for(j=0;j<nofabs;++j){
					outf2 << refseq[i].EnrPerTSS.peakbins[j] << '\t';
					outf2 << refseq[i].EnrPerTSS.nofpeakbins[j] << '\t';
					outf2 << refseq[i].EnrPerTSS.nonpeakbins[j] << '\t';
					outf2 << refseq[i].EnrPerTSS.nofnonpeakbins[j] << '\t';			
				}
				outf2 << endl;
			}
		}
	}



}
bool PromoterClass::AnnotateWithPromoters(string ichr, int ipos){ // Eliminate or include Promoter-promoter interactions
int i=0,j=0,k=0,z=0,geneindex=0;
int promid=-1, whichisoform=0;
bool pann=0,onprobe=1;

ichr.erase(ichr.find_last_not_of(" \n\r\t")+1); //trim the string
for(k=0;k<ChrNames_refseq.size();++k){
	if(ChrNames_refseq[k].compare(ichr.c_str())==0){ // Find the right index
		break;
	}
}
if(k==ChrNames_refseq.size())
	return 0;

for(z=0;z<refseq_indexes[k].size();++z){ //Iterate over all refseq genes on that chromosome
	for(j=0;j<refseq[refseq_indexes[k][z]].isoformpromotercoords.size();j++){ 
		if(abs(refseq[refseq_indexes[k][z]].isoformpromotercoords[j]-ipos)<=AssociateInteractions){ // If it is close enough to a promoter
			return 1; // eliminate promoter promoter inteactions
		}
	}
}
return 0;

}
void PromoterClass::AssociateProbeswithPromoters(ProbeSet& mm9prs){

int i=0,j=0,k=0,x=0;
int diff;
long int startsearch,endsearch;


for(i=0;i<NofPromoters;++i){
	startsearch=0;endsearch=-1;
	for(j=0;j<mm9prs.ChrNames.size();++j){
		if(refseq[i].chr==mm9prs.ChrNames[j]){
			startsearch=mm9prs.ChrRowStartIndexes[j];
			endsearch=mm9prs.ChrRowEndIndexes[j];
			break;
		}
	}
	x=0;
	for(j=startsearch;j<=endsearch;++j){
		for(k=0;k<refseq[i].isoformpromotercoords.size();++k){
			diff=mm9prs.MM9Probes[j].center-refseq[i].isoformpromotercoords[k];
			if(abs(diff)<AssociateInteractions){
				if(mm9prs.MM9Probes[j].annotated==0){
					mm9prs.MM9Probes[j].annotated=1; // Associate to a Promoter
					refseq[i].ProbeIDs.push_back(j);
				}
				else
					mm9prs.MM9Probes[j].conflicting_annotations=1;
			}
		}
		if(diff>BinSize)
			break;
	}
	if(i%10000==0)
		cout << i << "    Probes Associated with Promoters" << endl;
}

}

void PromoterClass::ReadBinCoverage(string InputFileNameBase,int ExperimentNum, RESitesClass& dpnIIsites, MappabilityClass& mappability){ // If the SAMFILE already processed, read the number of reads per bin
  //int i,j,k,t;
	string infilename,temp;
#ifdef UNIX
	infilename.append("/bubo/home/h20/pelin/3Cproj/bin/Detect_Interactions/inputFiles/");
#endif
	/*		
	infilename.append(InputFileNameBase);
	infilename.append("CoveragePerBin_proms.txt");
	ifstream infile(infilename.c_str());	
     
	for(int i=0;i<NofPromoters;++i){
		for(int k=0;k<refseq[i].isoformpromotercoords.size();++k){
			for(int t=0;t<5;t++)
				infile >> temp;
			for(int j=1;j<NumberofBins;++j){
			  int bincoordst = 0, recount = 0;
			  double mapp;
			  infile >> refseq[i].AllExperiments_IntBins[ExperimentNum].interactorbins[k][j];
			  //get the bin coordinate
			  if (refseq[i].strand == "+" ){
			    if (j < NumberofBins/2 )
			      bincoordst = refseq[i].isoformpromotercoords[k] - (NofInteractorBins-j)*BinSize; 
			    else
			      bincoordst = refseq[i].isoformpromotercoords[k] + (j-NofInteractorBins)*BinSize; 
			  }
			  else{
			    if ( j < NumberofBins/2)
			      bincoordst = refseq[i].isoformpromotercoords[k] + (NofInteractorBins-j)*BinSize; 
			    else
			      bincoordst = refseq[i].isoformpromotercoords[k] - (j- NofInteractorBins)*BinSize; 
			  }
			  if (refseq[i].AllExperiments_IntBins[ExperimentNum].interactorbins[k][j] > MinNumberofPairs){
			    mapp = mappability.GetMappability(refseq[i].chr, bincoordst);
			    if (mapp < 0.80)
			      refseq[i].AllExperiments_IntBins[ExperimentNum].interactorbins[k][j] = 0;
			    else{
			      recount = dpnIIsites.GetRESitesCount(refseq[i].chr, bincoordst);
			      if(recount != 0)
				refseq[i].AllExperiments_IntBins[ExperimentNum].interactorbins[k][j]/=recount;
			    }
			  }			 
			}
		}
		if (i%2000 == 0)
			cout << i << "   Promoter Coverage Read " << endl;
	}
	infile.close();
	*/	
	infilename.append(InputFileNameBase);
	infilename.append("NormalisedReadCounts_READTHIS.txt");
	ifstream infile(infilename.c_str());
	for(int i=0;i<NofPromoters;++i){
		for(int k=0;k<refseq[i].isoformpromotercoords.size();++k){
			for(int j=1;j<NumberofBins;++j){
				infile >> refseq[i].AllExperiments_IntBins[ExperimentNum].interactorbins[k][j];
			}
		}
	}
              
       
	
}
void PromoterClass::WriteBinCoverage(string OutputFileNameBase,int ExperimentNum){
	int i,j,k;
	string outfilename;
	
	outfilename.append("/bubo/proj/b2011029/bin/Detect_Interactions/outputFiles/");
	outfilename.append(OutputFileNameBase);
	outfilename.append("CoveragePerBin_proms.txt");
	ofstream outfile(outfilename.c_str());
	cout << refseq[3].RefSeqName << '\t' << refseq[3].chr << '\t' << refseq[3].isoformpromotercoords.size() << endl;
	for(i=0;i<NofPromoters;++i){
		for(k=0;k<refseq[i].isoformpromotercoords.size();++k){
			outfile << refseq[i].RefSeqName << '\t' << refseq[i].TranscriptName[k] << '\t' << refseq[i].chr << ":" << refseq[i].isoformpromotercoords[k]-5000 << "-" << refseq[i].isoformpromotercoords[k]+5000 << '\t';
			outfile << refseq[i].strand << '\t' << refseq[i].expression[CellType] << '\t';
			for(j=1;j<NumberofBins;++j)
				outfile << refseq[i].AllExperiments_IntBins[ExperimentNum].interactorbins[k][j] << '\t';
			outfile << endl;
		}
	}
	outfile.close();
}
