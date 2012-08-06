struct IntStruct{
vector <int> NofInteractionsAssociatedwithPeaks; //[0] = peak+, int+ [1] = peak+,int-, [2]=peak-,int+, [3]=peak-,int-
vector < vector <int> > distances; // Distance of the above classes
};

class DetectEnrichedBins{
public:
	IntStruct *TSS_Ints; //Int: interactions
	IntStruct *NG_Ints; 
	IntStruct *TSS_Ints_meta; // all peak data combined
	IntStruct *NG_Ints_meta;
	int *EnrichedIntBins; // How many enriched bins per positions
	void InitialiseData();
	double partialfactorial(double,double);
	void CalculateEnrichments(IntStruct,ofstream&);
	void DetectEnrichedInteractionBins(PromoterClass&,double*,double*,string,int); // Find bins that have above background level interactions 
	void DetectEnrichedInteractionBins_NegCtrls(NegCtrlClass&,double*,double*,string,int); // Find bins that have above background level interactions 
	void AssociatePeaksWithIntBins(PromoterClass&,string,string,int,int,int);
	void AssociatePeaksWithIntBins_NegCtrls(NegCtrlClass&,string,string,int,int);
	void PrintMetaAssociationwithPeaks(PromoterClass&,NegCtrlClass&,string,int);

	~DetectEnrichedBins();
};

DetectEnrichedBins::~DetectEnrichedBins(){

}
void DetectEnrichedBins::InitialiseData(void){
EnrichedIntBins = new int[NumberofBins];

}
double DetectEnrichedBins::partialfactorial(double n, double r){

	unsigned int i=0;
	double factorial;

	factorial=1;
	if(r==0){
		for(i=n;i>0;--i)
			factorial*=i;

		return (factorial);
	}
	else{
		for(i=n;i>r;--i)
		factorial*=i;
	}

	return (factorial);

}
void DetectEnrichedBins::DetectEnrichedInteractionBins_NegCtrls(NegCtrlClass &ngs,double* BGfreqs,double *stdevs,string BaseFileName,int whichexperiment){
	for(unsigned int i=0;i<NumberofBins;++i)
		EnrichedIntBins[i]=0;

	string FileName;
	FileName.append(BaseFileName);
	FileName.append("_NofINTPerBin_negctrls_");
	FileName.append(".txt");
	ofstream outf1(FileName.c_str());

for(unsigned int i=0;i<ngs.NofNegCtrls;++i){
			for(unsigned int j=1;j<NumberofBins;++j){
				double z_score,q,p_value;
				z_score=((ngs.negctrls[i].AllExperiments_IntBins[whichexperiment].interactorbins[0][j]-BGfreqs[j])/stdevs[j]);
				q=alglib::normaldistribution(z_score);
				p_value=1-q;
				if(p_value<SignificanceThreshold && ngs.negctrls[i].AllExperiments_IntBins[whichexperiment].interactorbins[0][j]>MinNumberofPairs){
					++EnrichedIntBins[j];
					ngs.negctrls[i].AllExperiments_IntBins[whichexperiment].Interactor[0][j]=p_value;
					ngs.negctrls[i].AllExperiments_IntBins[whichexperiment].EnrichedBins[0].push_back(j);				
				}
				else
					ngs.negctrls[i].AllExperiments_IntBins[whichexperiment].Interactor[0][j]=1;
			}
}
	for(int j=1;j<NumberofBins;++j)
		outf1 << EnrichedIntBins[j] << '\t';
	outf1 << endl;

}
void DetectEnrichedBins::DetectEnrichedInteractionBins(PromoterClass &prs,double* BGfreqs,double *stdevs,string BaseFileName,int whichexperiment){
	for(int i=0;i<NumberofBins;++i)
		EnrichedIntBins[i]=0;

	string FileName;
	FileName.append(BaseFileName);
	FileName.append("_NofINTPerBin");
	FileName.append(".txt");
	ofstream outf1(FileName.c_str());
//#pragma omp parallel for default(shared)
for(int i=0;i<prs.NofPromoters;++i){
	for(int k=0;k<prs.refseq[i].isoformpromotercoords.size();++k){
		for(int j=1;j<NumberofBins;++j){
			double z_score,q,p_value;
			z_score=((prs.refseq[i].AllExperiments_IntBins[whichexperiment].interactorbins[k][j]-BGfreqs[j])/stdevs[j]);
			q=alglib::normaldistribution(z_score);
			p_value=1-q;
			if(p_value<SignificanceThreshold && prs.refseq[i].AllExperiments_IntBins[whichexperiment].interactorbins[k][j]>MinNumberofPairs){
				++EnrichedIntBins[j];
				prs.refseq[i].AllExperiments_IntBins[whichexperiment].Interactor[k][j]=p_value;
				prs.refseq[i].AllExperiments_IntBins[whichexperiment].EnrichedBins[k].push_back(j);				
			}
			else
				prs.refseq[i].AllExperiments_IntBins[whichexperiment].Interactor[k][j]=1;
		}
	}

}
for(int j=1;j<NumberofBins;++j)
	outf1 << EnrichedIntBins[j] << '\t';
outf1 << endl;
}
void DetectEnrichedBins::AssociatePeaksWithIntBins(PromoterClass& prs, string BaseFileName,string abname,int CellType,int abindex,int ExperimentIndex){
	int i,j,k,l,indx=0;
	string FileName;
	FileName.append(BaseFileName);
	FileName.append(abname);
	FileName.append(".txt");
	ofstream outf1(FileName.c_str());

TSS_Ints = new IntStruct[NumberofGenes];

outf1 << "Gene Name" << '\t' << "Transcript_Name" << '\t' << "Expression" << '\t' << "UCSC" << '\t' << "Peak w Interaction" << '\t' << "Peak w/o Interaction" << '\t' << "No peak but Interaction" << '\t' << "No peak no interaction" << '\t';
outf1 << "Number of Peak Bins with Interaction/Number of Bins with Interaction" << '\t' << "Number of Peak Bins/Total Number of Bins" << '\t' << "Number of Bins with Interactions" << '\t';
outf1 << "Probability of occurence by chance" << '\t' << "Enrichment" << endl;
for(i=0;i<prs.NofPromoters;++i){
	for(l=0;l<prs.refseq[i].isoformpromotercoords.size();++l){
		TSS_Ints[indx].NofInteractionsAssociatedwithPeaks.resize(4);
			for(j=1;j<NumberofBins;++j){
				if(prs.refseq[i].AllPeaks_PeakBins[abindex].PeakBins[l][j]>MinNumberofPeaksPerBin){
					if(prs.refseq[i].AllExperiments_IntBins[ExperimentIndex].Interactor[l][j]<SignificanceThreshold)
						++TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[0];
					else
						++TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[1];
				}
				else{
					if(prs.refseq[i].AllExperiments_IntBins[ExperimentIndex].Interactor[l][j]<SignificanceThreshold)
						++TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[2];
					else
						++TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[3];
				}
			}
			outf1 << prs.refseq[i].RefSeqName << '\t' << prs.refseq[i].TranscriptName[l] << '\t' << prs.refseq[i].expression[CellType] << '\t';
			outf1 << prs.refseq[i].chr << ":" << prs.refseq[i].isoformpromotercoords[l]-5000 << "-" << prs.refseq[i].isoformpromotercoords[l]+5000 << '\t';
			for(k=0;k<4;++k)
				outf1 << TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[k] << '\t';
			CalculateEnrichments(TSS_Ints[indx],outf1);
			++indx;
	}
}
	delete[] TSS_Ints;
}
void DetectEnrichedBins::AssociatePeaksWithIntBins_NegCtrls(NegCtrlClass& ngs,string BaseFileName,string abname,int abindex,int ExperimentIndex){
	int i,j,k,indx=0;
	string FileName;
	FileName.append(BaseFileName);
	FileName.append(abname);
	FileName.append("_negctrls");
	FileName.append(".txt");
	ofstream outf1(FileName.c_str());

	NG_Ints = new IntStruct[ngs.NofNegCtrls];

outf1 << "Negative Control ID" << '\t' << "Peak w Interaction" << '\t' << "Peak w/o Interaction" << '\t' << "No peak but Interaction" << '\t' << "No peak no interaction" << endl;
outf1 << "Number of Peak Bins with Interaction/Number of Bins with Interaction" << '\t' << "Number of Peak Bins/Total Number of Bins" << '\t' << "Number of Bins with Interactions" << '\t';
outf1 << "Probability of occurence by chance" << '\t' << "Enrichment" << endl;

for(i=0;i<ngs.NofNegCtrls;++i){
	NG_Ints[indx].NofInteractionsAssociatedwithPeaks.resize(4);
			for(j=1;j<NumberofBins;++j){
				if(ngs.negctrls[i].AllPeaks_PeakBins[abindex].PeakBins[0][j]>0){
					if(ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].Interactor[0][j]<SignificanceThreshold)
						++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[0];
					else
						++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[1];
				}
				else{
					if(ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].Interactor[0][j]<SignificanceThreshold)
						++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[2];
					else	
						++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[3];
				}
			}
			outf1 << i  << '\t';
			for(k=0;k<4;++k)
				outf1 << NG_Ints[indx].NofInteractionsAssociatedwithPeaks[k] << '\t';
			CalculateEnrichments(NG_Ints[indx],outf1);
			++indx;
}
	delete[] NG_Ints;


}
void DetectEnrichedBins::CalculateEnrichments(IntStruct interactions,ofstream& outf){
	double NumberofInteractions=double((interactions.NofInteractionsAssociatedwithPeaks[0]+interactions.NofInteractionsAssociatedwithPeaks[2]));
	double a;
	//a = Number of Peak Bins with Interaction/Number of Bins with Interaction
	if(NumberofInteractions>0)
		a = double(interactions.NofInteractionsAssociatedwithPeaks[0])/(NumberofInteractions);
	else
		a = 0;
	// b = Number of Peak Bins/Total Number of Bins
	double b = double((interactions.NofInteractionsAssociatedwithPeaks[0]+interactions.NofInteractionsAssociatedwithPeaks[1]))/(NumberofBins-1);
	//c = Number of Bins with Interactions
	double c = double((interactions.NofInteractionsAssociatedwithPeaks[0]+interactions.NofInteractionsAssociatedwithPeaks[2]));
	// p = Probability of occurence by chance
	double t1 = (pow(b,(double (interactions.NofInteractionsAssociatedwithPeaks[0])))*pow((1-b),(double (interactions.NofInteractionsAssociatedwithPeaks[2]))));
	double r1 = (double(c-interactions.NofInteractionsAssociatedwithPeaks[0]));
	double t2 = partialfactorial(c,r1);
	double t3=(partialfactorial((double(interactions.NofInteractionsAssociatedwithPeaks[0])),0));
	double p = t1*(t2/t3);
	//			(POWER(I2,D2)*POWER((1-I2),F2))*(FACT(J2)/(FACT(D2)*(FACT(F2))))
	//			(POWER(b,[0])*POWER((1-b),[2]))*(FACT(c)/(FACT([0])*(FACT([3]))))
	// e = Enrichment
	double e = a/p;
	outf << a << '\t' << b << '\t' << c << '\t' << p << '\t' << e << endl;

}
void DetectEnrichedBins::PrintMetaAssociationwithPeaks(PromoterClass& prs,NegCtrlClass& ngs,string BaseFileName,int ExperimentIndex){

	int i,j,k,l,indx=0,ab;

TSS_Ints = new IntStruct[NumberofGenes];
string filen1;
bool pi_bin,p_bin;
filen1.append(BaseFileName);
filen1.append("MetaAssociationwithPeaks_Promoters.txt");
ofstream outf1(filen1.c_str());

outf1 << "Gene Name" << '\t' << "Transcript_Name" << '\t' << "Expression" << '\t' << "Number of Probes" << '\t' << "Peak w Interaction" << '\t' << "Peak w/o Interaction" << '\t' << "No peak but Interaction" << '\t' << "No peak no interaction" << '\t';
outf1 << "Number of Peak Bins with Interaction/Number of Bins with Interaction" << '\t' << "Number of Peak Bins/Total Number of Bins" << '\t' << "Number of Bins with Interactions" << '\t';
outf1 << "Probability of occurence by chance" << '\t' << "Enrichment" << endl;

for(i=0;i<prs.NofPromoters;++i){
	for(l=0;l<prs.refseq[i].isoformpromotercoords.size();++l){
		TSS_Ints[indx].NofInteractionsAssociatedwithPeaks.resize(4);
		for(j=0;j<4;++j)
			TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[j]=0;
		for(j=1;j<NumberofBins;++j){
			pi_bin=0; p_bin=0;
			if(prs.refseq[i].AllExperiments_IntBins[ExperimentIndex].Interactor[l][j]<SignificanceThreshold){
				for(ab=0;ab<NumberofPeakFiles;++ab){
					if(prs.refseq[i].AllPeaks_PeakBins[ab].PeakBins[l][j]>0){
						TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[0]++; //Peak with Interaction
						pi_bin=1;
						break;
					}
				}
				if(!pi_bin)
					TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[2]++;
			}
			else{
				for(ab=0;ab<NumberofPeakFiles;++ab){
					if(prs.refseq[i].AllPeaks_PeakBins[ab].PeakBins[l][j]>0){
						TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[1]++; //Peak with Interaction
						p_bin=1;
						break;
					}
				}
				if(!p_bin)
					TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[3]++;
			}

		}
		outf1 << prs.refseq[i].RefSeqName << '\t' << prs.refseq[i].TranscriptName[l] << '\t' << prs.refseq[i].expression[CellType] << '\t' << prs.refseq[i].ProbeIDs.size() << '\t';
		for(k=0;k<4;++k)
			outf1 << TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[k] << '\t';
		CalculateEnrichments(TSS_Ints[indx],outf1);

		++indx;
	}
}
delete[] TSS_Ints;
outf1.close();

NG_Ints = new IntStruct[ngs.NofNegCtrls];
string filen2;
filen2.append(BaseFileName);
filen2.append("MetaAssociationwithPeaks_NegCtrls.txt");

ofstream outf2(filen2.c_str());
indx=0;
outf2 << "Negative Control ID" << '\t' << "Peak w Interaction" << '\t' << "Peak w/o Interaction" << '\t' << "No peak but Interaction" << '\t' << "No peak no interaction" << endl;
outf2 << "Number of Peak Bins with Interaction/Number of Bins with Interaction" << '\t' << "Number of Peak Bins/Total Number of Bins" << '\t' << "Number of Bins with Interactions" << '\t';
outf2 << "Probability of occurence by chance" << '\t' << "Enrichment" << endl;

for(i=0;i<ngs.NofNegCtrls;++i){
	NG_Ints[indx].NofInteractionsAssociatedwithPeaks.resize(4);
	for(j=0;j<4;++j)
		NG_Ints[indx].NofInteractionsAssociatedwithPeaks[j]=0;
	for(j=1;j<NumberofBins;++j){
		pi_bin=0; p_bin=0;
		if(ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].Interactor[0][j]<SignificanceThreshold){
			for(ab=0;ab<NumberofPeakFiles;++ab){
					if(ngs.negctrls[i].AllPeaks_PeakBins[ab].PeakBins[0][j]>0){
						if(ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].Interactor[0][j]<SignificanceThreshold){
							++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[0];
							pi_bin=1;
							break;
						}
					}
			}
			if(!pi_bin)
				++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[2];
		}
		else{
			for(ab=0;ab<NumberofPeakFiles;++ab){
				if(ngs.negctrls[i].AllPeaks_PeakBins[ab].PeakBins[0][j]>0){
					if(ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].Interactor[0][j]<SignificanceThreshold){
						++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[1];
						p_bin=1;
						break;
					}
				}
			}
			if(!p_bin)
				++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[3];
		}
	}
		outf2 << i  << '\t';
		for(k=0;k<4;++k)
			outf2 << NG_Ints[indx].NofInteractionsAssociatedwithPeaks[k] << '\t';
		CalculateEnrichments(NG_Ints[indx],outf2);
		++indx;
}
delete[] NG_Ints;
outf2.close();
}