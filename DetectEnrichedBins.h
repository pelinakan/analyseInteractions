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
	vector < vector < int > > ClusteredEnrichedBins; // Enriched bins are clustered to minimise false positives
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
	FileName.append("/bubo/proj/b2011029/bin/Detect_Interactions/outputFiles/");
	FileName.append(BaseFileName);
	FileName.append("_NofINTPerBin_negctrls_");
	FileName.append(".txt");
	ofstream outf1(FileName.c_str());

	vector < int > TempEnrichedBins;

for(unsigned int i=0;i<ngs.NofNegCtrls;++i){
			for(unsigned int j=1;j<NumberofBins;++j){
		  if(j < ((NumberofBins/2)-CloseBins) || j > ((NumberofBins/2)+CloseBins)){ // do not look for the neighboring 10 kb upstream and downstream
					double z_score,q,p_value;
					z_score=((ngs.negctrls[i].AllExperiments_IntBins[whichexperiment].interactorbins[0][j]-BGfreqs[j])/stdevs[j]);
					q=alglib::normaldistribution(z_score);
					p_value=1-q;
					if(p_value<SignificanceThreshold && ngs.negctrls[i].AllExperiments_IntBins[whichexperiment].interactorbins[0][j]>MinNumberofPairs){
						++EnrichedIntBins[j];
						ngs.negctrls[i].AllExperiments_IntBins[whichexperiment].Interactor[0][j]=p_value;
						TempEnrichedBins.push_back(j);				
					}
					else
						ngs.negctrls[i].AllExperiments_IntBins[whichexperiment].Interactor[0][j]=1;
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
	
//	GCContent GC;
//	GC.InitialiseVars();
//	double gc1, gc2;

	string FileName,FileName2,FileName3, FileName4;
	FileName.append("/bubo/proj/b2011029/bin/Detect_Interactions/outputFiles/");
	FileName.append(BaseFileName);
	FileName.append("_NofINTPerBin_");
	FileName.append(FilterInteractions);
	FileName.append(".txt");
	ofstream outf1(FileName.c_str());

	FileName3.append("/bubo/proj/b2011029/bin/Detect_Interactions/outputFiles/");
	FileName3.append(BaseFileName);
	FileName3.append("_NormalisedReadCounts_");
	FileName.append(FilterInteractions);
	FileName3.append(".txt");
	ofstream outf3(FileName3.c_str());
/*
	FileName4.append(BaseFileName);
	FileName4.append("_GCcontent_PE");
	FileName4.append(".txt");
	ofstream outf4(FileName4.c_str());
*/
	vector < int > TempEnrichedBins;

for(int i=0;i<prs.NofPromoters;++i){
	for(int k=0;k<prs.refseq[i].isoformpromotercoords.size();++k){
		for(int j=1;j<NumberofBins;++j){
		  if(j < ((NumberofBins/2)-CloseBins) || j > ((NumberofBins/2)+CloseBins)){ // do not look for the neighboring 10 kb upstream and downstream
				double z_score,q,p_value;
				z_score=((prs.refseq[i].AllExperiments_IntBins[whichexperiment].interactorbins[k][j]-BGfreqs[j])/stdevs[j]);
				q=alglib::normaldistribution(z_score);
				p_value=1-q;
				if((p_value<SignificanceThreshold) && prs.refseq[i].AllExperiments_IntBins[whichexperiment].interactorbins[k][j] >= MinNumberofPairs){
					int bincoordst = 0;
					if (prs.refseq[i].strand == "+" ){
						if (j < NumberofBins/2 )
							bincoordst = prs.refseq[i].isoformpromotercoords[k] - (NofInteractorBins-j)*BinSize; 
						else
							bincoordst = prs.refseq[i].isoformpromotercoords[k] + (j-NofInteractorBins)*BinSize; 
					}
					else{
						if ( j < NumberofBins/2)
							bincoordst = prs.refseq[i].isoformpromotercoords[k] + (NofInteractorBins-j)*BinSize; 
						else
							bincoordst = prs.refseq[i].isoformpromotercoords[k] - (j- NofInteractorBins)*BinSize; 
					}
					bool pp = prs.AnnotateWithPromoters(prs.refseq[i].chr, bincoordst); //Promoter-Promoter Int ?
					if (pp){
//						gc1 = GC.GetGCContent(prs.refseq[i].chr,prs.refseq[i].isoformpromotercoords[k]);
//						gc2 = GC.GetGCContent(prs.refseq[i].chr,bincoordst);
//						outf4 << gc1 << '\t' << gc2 << '\t' << prs.refseq[i].AllExperiments_IntBins[whichexperiment].interactorbins[k][j] << endl;
						prs.refseq[i].AllExperiments_IntBins[whichexperiment].Interactor[k][j] = p_value;
						TempEnrichedBins.push_back(j);
					}
					else
						prs.refseq[i].AllExperiments_IntBins[whichexperiment].Interactor[k][j] = 1;
				}
				else
					prs.refseq[i].AllExperiments_IntBins[whichexperiment].Interactor[k][j] = 1;
			}
			else
				prs.refseq[i].AllExperiments_IntBins[whichexperiment].Interactor[k][j] = 1;

			outf3 << prs.refseq[i].AllExperiments_IntBins[whichexperiment].interactorbins[k][j] <<  '\t';
		}
		outf3 << endl;
	}
}
for(int j=1;j<NumberofBins;++j)
	outf1 << EnrichedIntBins[j] << '\t';
outf1 << endl;
}
void DetectEnrichedBins::AssociatePeaksWithIntBins(PromoterClass& prs, string BaseFileName,string abname,int CellType,int abindex,int ExperimentIndex){
	int i,j,k,l,indx=0;
	string FileName,FileName2,FileName3;

	FileName.append("/bubo/home/h20/pelin/3Cproj/bin/Detect_Interactions/outputFiles/");
	FileName.append(BaseFileName);
	FileName.append(abname);
	FileName.append("_");
	FileName.append(FilterInteractions);
	FileName.append(".txt");
	ofstream outf1(FileName.c_str());

	FileName2.append("/bubo/home/h20/pelin/3Cproj/bin/Detect_Interactions/outputFiles/");
	FileName2.append(BaseFileName);
	FileName2.append(abname);
	FileName2.append("_");
	FileName2.append(FilterInteractions);
	FileName2.append("_IntScore_PeakOverlap.txt");
	ofstream outf2(FileName2.c_str());

	FileName3.append("/bubo/home/h20/pelin/3Cproj/bin/Detect_Interactions/outputFiles/");
	FileName3.append(BaseFileName);
	FileName3.append(abname);
	FileName3.append("_");
	FileName3.append(FilterInteractions);
	FileName3.append("_IntScore_NoPeakOverlap.txt");
	ofstream outf3(FileName3.c_str());

	//GCContent GC;
	//	GC.InitialiseVars();
	double gc1, gc2;


TSS_Ints = new IntStruct[NumberofGenes];
double sum1 = 0, sum2 = 0, sum3 = 0;
int bindex = 0;
int bincoordst = 0;


outf1 << "Gene Name" << '\t' << "Transcript_Name" << '\t' << "Expression" << '\t' << "UCSC" << '\t' << "Peak w Interaction" << '\t' << "Peak w/o Interaction" << '\t' << "No peak but Interaction" << '\t' << "No peak no interaction" << '\t' << "No peak no interaction" << '\t';
outf1 << "Number of Peak Bins with Interaction/Number of Bins with Interaction" << '\t' << "Number of Peak Bins/Total Number of Bins" << '\t' << "Number of Bins with Interactions" << '\t';
outf1 << "Probability of occurence by chance" << '\t' << "Enrichment" << endl;
for(i=0;i<prs.NofPromoters;++i){
	for(l = 0; l < prs.refseq[i].isoformpromotercoords.size();++l){
		TSS_Ints[indx].NofInteractionsAssociatedwithPeaks.resize(4);
		for(j=0;j<4;++j)
			TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[j]=0;
		for(j = 1; j < NumberofBins;++j){
			bool peakoverlap = 0;
			if(prs.refseq[i].AllExperiments_IntBins[ExperimentIndex].Interactor[l][j] < SignificanceThreshold){
				if(prs.refseq[i].AllPeaks_PeakBins[abindex].PeakBins[l][j] > MinNumberofPeaksPerBin){
					++TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[0];
					outf2 <<  prs.refseq[i].AllExperiments_IntBins[ExperimentIndex].interactorbins[l][j] << endl;
				}
				else{
					++TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[2];
					outf3 << prs.refseq[i].AllExperiments_IntBins[ExperimentIndex].interactorbins[l][j] << endl;
				}
			}
			else{
				if(prs.refseq[i].AllPeaks_PeakBins[abindex].PeakBins[l][j] > MinNumberofPeaksPerBin)
					++TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[1];
				else
					++TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[3];
			}
		}
		outf1 << prs.refseq[i].RefSeqName << '\t' << prs.refseq[i].TranscriptName[l] << '\t' << prs.refseq[i].expression[CellType] << '\t';
		outf1 << prs.refseq[i].chr << ":" << prs.refseq[i].isoformpromotercoords[l]-5000 << "-" << prs.refseq[i].isoformpromotercoords[l]+5000 << '\t';
		for(k=0;k<4;++k)
			outf1 << TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[k] << '\t';
		CalculateEnrichments(TSS_Ints[indx],outf1);
		sum1 += TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[0];
		sum2 += TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[1];
		sum3 += TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[2];
		++indx;
	}
}
po_file << abname << '\t' << sum1 << '\t' << sum2 << '\t' << sum3 << '\t' <<  sum1+sum3 << '\t' << sum1/(sum1+sum3) << endl;
	delete[] TSS_Ints;

}
void DetectEnrichedBins::AssociatePeaksWithIntBins_NegCtrls(NegCtrlClass& ngs,string BaseFileName,string abname,int abindex,int ExperimentIndex){
	int i,j,k,indx=0;
	string FileName;

	FileName.append("/bubo/proj/b2011029/bin/Detect_Interactions/outputFiles/");
	FileName.append(BaseFileName);
	FileName.append(abname);
	FileName.append("_negctrls");
	FileName.append(".txt");
	ofstream outf1(FileName.c_str());

	NG_Ints = new IntStruct[ngs.NofNegCtrls];

	double sum1 = 0, sum2 = 0, sum3 = 0;
outf1 << "Negative Control ID" << '\t' << "Peak w Interaction" << '\t' << "Peak w/o Interaction" << '\t' << "No peak but Interaction" << '\t' << "No peak no interaction" << '\t';
outf1 << "Number of Peak Bins with Interaction/Number of Bins with Interaction" << '\t' << "Number of Peak Bins/Total Number of Bins" << '\t' << "Number of Bins with Interactions" << '\t';
outf1 << "Probability of occurence by chance" << '\t' << "Enrichment" << endl;

for(i=0;i<ngs.NofNegCtrls;++i){
	NG_Ints[indx].NofInteractionsAssociatedwithPeaks.resize(4);
			for(j=1;j<NumberofBins;++j){
				if(ngs.negctrls[i].AllPeaks_PeakBins[abindex].PeakBins[0][j] > MinNumberofPeaksPerBin){
					if(ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].Interactor[0][j]<SignificanceThreshold && ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].interactorbins[0][j] >= MinNumberofPairs)
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
			sum1 += NG_Ints[indx].NofInteractionsAssociatedwithPeaks[0];
			sum2 += NG_Ints[indx].NofInteractionsAssociatedwithPeaks[1];
			sum3 += NG_Ints[indx].NofInteractionsAssociatedwithPeaks[2];
			++indx;
}
	delete[] NG_Ints;
	
	po_file << abname << "_NegCtrls" << '\t' << sum1 << '\t' << sum2 << '\t' << sum3 << '\t' <<  sum1+sum3 << '\t' << sum1/(sum1+sum3) << endl;

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
	double sum1 = 0, sum2 = 0, sum3 = 0;

TSS_Ints = new IntStruct[NumberofGenes];
string filen1;

filen1.append("/bubo/proj/b2011029/bin/Detect_Interactions/outputFiles/");
filen1.append(BaseFileName);
 filen1.append(FilterInteractions);
 filen1.append("_");
filen1.append("MetaAssociationwithPeaks_Promoters.txt");
ofstream outf1(filen1.c_str());

outf1 << "Gene Name" << '\t' << "Transcript_Name" << '\t' << "Expression" << '\t' << "Number of Probes" << '\t' << "Peak w Interaction" << '\t' << "Peak w/o Interaction" << '\t' << "No peak but Interaction" << '\t';
outf1 << "Number of Peak Bins with Interaction/Number of Bins with Interaction" << '\t' << "Number of Peak Bins/Total Number of Bins" << '\t' << "Number of Bins with Interactions" << '\t';
outf1 << "Probability of occurence by chance" << '\t' << "Enrichment" << endl;

string FileName2, FileName3;
        FileName2.append("/bubo/proj/b2011029/bin/Detect_Interactions/outputFiles/");
	FileName2.append(BaseFileName);
	FileName2.append("Meta_IntScore_PeakOverlap.txt");
	ofstream outfP(FileName2.c_str());

	FileName3.append("/bubo/proj/b2011029/bin/Detect_Interactions/outputFiles/");
	FileName3.append(BaseFileName);
	FileName3.append("Meta_IntScore_NoPeakOverlap.txt");
	ofstream outfNP(FileName3.c_str());

for(i=0;i<prs.NofPromoters;++i){
	for(l=0;l<prs.refseq[i].isoformpromotercoords.size();++l){
		TSS_Ints[indx].NofInteractionsAssociatedwithPeaks.resize(4);
		for(j=0;j<4;++j)
			TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[j]=0;
		for(j = 1; j < NumberofBins;++j){
			bool peakbin = 0;
			for(ab=0;ab<NumberofPeakFiles;++ab){
				if(prs.refseq[i].AllPeaks_PeakBins[ab].PeakBins[l][j] > MinNumberofPeaksPerBin){
					peakbin = 1;
					break;
				}
			}
			if(prs.refseq[i].AllExperiments_IntBins[ExperimentIndex].Interactor[l][j] < SignificanceThreshold){
				if (peakbin){
					++TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[0];
					outfP <<  prs.refseq[i].AllExperiments_IntBins[ExperimentIndex].interactorbins[l][j] << endl;
				}
				else{
					++TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[2];
					outfNP <<  prs.refseq[i].AllExperiments_IntBins[ExperimentIndex].interactorbins[l][j] << endl;
				}
			}
			else{
				if (peakbin)
					++TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[1];
				else
					++TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[3];

			}
		}
		sum1 += TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[0];
		sum2 += TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[1];
		sum3 += TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[2];
		
		outf1 << prs.refseq[i].RefSeqName << '\t' << prs.refseq[i].TranscriptName[l] << '\t' << prs.refseq[i].expression[CellType] << '\t' << prs.refseq[i].ProbeIDs.size() << '\t';
		for(k=0;k<4;++k)
			outf1 << TSS_Ints[indx].NofInteractionsAssociatedwithPeaks[k] << '\t';
		CalculateEnrichments(TSS_Ints[indx],outf1);
		++indx;
	}
}

delete[] TSS_Ints;
outf1.close();
cout << "Meta Overlap of Peaks with Interactions   " << (sum1/(sum1+sum3)) << endl;
cout << "Total Number of Interactions    " << sum1+sum3 << endl;
po_file << "Meta" <<  '\t' << sum1 <<  '\t' << sum2 <<  '\t' << sum3  << '\t' << sum1+sum3 << '\t' << (sum1/(sum1+sum3)) << endl;

NG_Ints = new IntStruct[ngs.NofNegCtrls];
string filen2;

filen2.append("/bubo/proj/b2011029/bin/Detect_Interactions/outputFiles/");
filen2.append(BaseFileName);
filen2.append("MetaAssociationwithPeaks_NegCtrls.txt");

ofstream outf2(filen2.c_str());
indx=0;
outf2 << "Negative Control ID" << '\t' << "Peak w Interaction" << '\t' << "Peak w/o Interaction" << '\t' << "No peak but Interaction" << '\t' << "No peak no interaction" << '\t';
outf2 << "Number of Peak Bins with Interaction/Number of Bins with Interaction" << '\t' << "Number of Peak Bins/Total Number of Bins" << '\t' << "Number of Bins with Interactions" << '\t';
outf2 << "Probability of occurence by chance" << '\t' << "Enrichment" << endl;

sum1 = 0; sum2 = 0; sum3 = 0;
for(i=0;i<ngs.NofNegCtrls;++i){
		NG_Ints[indx].NofInteractionsAssociatedwithPeaks.resize(4);
		for(j=0;j<4;++j)
			NG_Ints[indx].NofInteractionsAssociatedwithPeaks[j]=0;
		for(j = 1; j < NumberofBins;++j){
			bool peakbin = 0;
			for(ab=0;ab<NumberofPeakFiles;++ab){
				if(ngs.negctrls[i].AllPeaks_PeakBins[ab].PeakBins[0][j] > MinNumberofPeaksPerBin){
					peakbin = 1;
					break;
				}
			}
			if(ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].Interactor[0][j] < SignificanceThreshold){
				if (peakbin){
					++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[0];
				}
				else{
					++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[2];
				}
			}
			else{
				if (peakbin)
					++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[1];
				else
					++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[3];

			}
		}
		sum1 += NG_Ints[indx].NofInteractionsAssociatedwithPeaks[0];
		sum2 += NG_Ints[indx].NofInteractionsAssociatedwithPeaks[1];
		sum3 += NG_Ints[indx].NofInteractionsAssociatedwithPeaks[2];
		
		outf1 <<  i  << '\t';
		for(k=0;k<4;++k)
			outf1 << NG_Ints[indx].NofInteractionsAssociatedwithPeaks[k] << '\t';
		CalculateEnrichments(NG_Ints[indx],outf1);
		++indx;

}
delete[] NG_Ints;
po_file << "Meta_NegCtrls" <<  '\t' << sum1 <<  '\t' << sum2 <<  '\t' << sum3  << '\t' << sum1+sum3 << '\t' << (sum1/(sum1+sum3)) << endl;

outf2.close();
}

/*
	NG_Ints[indx].NofInteractionsAssociatedwithPeaks.resize(4);
	for(j=0;j<4;++j)
		NG_Ints[indx].NofInteractionsAssociatedwithPeaks[j]=0;
	for(j=1;j<NumberofBins;++j){
		bool pi_bin=0, p_bin=0;
		if(ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].Interactor[0][j]<SignificanceThreshold && ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].interactorbins[0][j] >= MinNumberofPairs){
			for(ab=0;ab<NumberofPeakFiles;++ab){
				if(ngs.negctrls[i].AllPeaks_PeakBins[ab].PeakBins[0][j]> MinNumberofPeaksPerBin){
						if(ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].Interactor[0][j]<SignificanceThreshold && ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].interactorbins[0][j] >= MinNumberofPairs){
							++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[0];						   
							pi_bin=1;
							break;
						}
					}
			}
			if (pi_bin)
			  sum2++;
			else
				++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[2];
		}
		else{
			for(ab=0;ab<NumberofPeakFiles;++ab){
				if(ngs.negctrls[i].AllPeaks_PeakBins[ab].PeakBins[0][j]>0){
					if(ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].Interactor[0][j]<SignificanceThreshold && ngs.negctrls[i].AllExperiments_IntBins[ExperimentIndex].interactorbins[0][j] >= MinNumberofPairs){
						++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[1];
						p_bin=1;
						break;
					}
				}
			}
			if (pi_bin)
			  sum2++;
			else
				++NG_Ints[indx].NofInteractionsAssociatedwithPeaks[3];
		}
	}
		outf2 << i  << '\t';
		for(k=0;k<4;++k)
			outf2 << NG_Ints[indx].NofInteractionsAssociatedwithPeaks[k] << '\t';
		CalculateEnrichments(NG_Ints[indx],outf2);
		sum1 += NG_Ints[indx].NofInteractionsAssociatedwithPeaks[0];
		//sum2 += NG_Ints[indx].NofInteractionsAssociatedwithPeaks[1];
		sum3 += NG_Ints[indx].NofInteractionsAssociatedwithPeaks[2];

		++indx;
 */
