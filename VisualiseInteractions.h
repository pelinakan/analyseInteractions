struct ColorStruct{
	int color[3];
	string colorname;
};
class VisualiseInteractions{
ColorStruct *Colors;
public:
void WriteGFFFiles(PromoterClass&,ofstream&);
void WriteBEDFiles(PromoterClass&,ofstream&);

private:
	void WritePromoterLine(PromoterClass&,ofstream&,int,int,int);
	void WriteEnhancerLine(PromoterClass&,ofstream&,int,int,int,int);
	void InitialiseColorMatrix(void);
	int DetermineColor(PromoterClass&,int,int,int);	
	void WriteBEDLine(PromoterClass&,ofstream&,int,int,int,int);
};
void VisualiseInteractions::InitialiseColorMatrix(void){
#ifdef UNIX
	ifstream cfile("/bubo/proj/b2011029/bin/3CAnalysis/ColorCoding.txt");
#endif	
	int index;
	string temp;
#ifdef WINDOWS
	ifstream cfile("ColorCoding.txt");
#endif
	Colors = new ColorStruct[7];
	for(int i=0;i<7;++i){
		cfile >> temp >> index;
		cfile >> Colors[index].colorname >> Colors[index].color[0]	>> Colors[index].color[1] >> Colors[index].color[2];
	}

}
void VisualiseInteractions::WriteGFFFiles(PromoterClass& prs,ofstream& GFFFile){

int gffindex;
int ColorIndex = 1;
for(int i=0;i<prs.NofPromoters;++i){
	gffindex=0;
	//		ColorIndex=DetermineColor(prs,i,j,b);
	//		if(ColorIndex!=-1){ // If there is an interaction
	for(int j = 0; j < prs.refseq[i].AllExperiments_IntBins[0].ClusteredEnrichedBins.size();++j){
		if( prs.refseq[i].AllExperiments_IntBins[0].ClusteredEnrichedBins[j].size() >= 1){
			WritePromoterLine(prs,GFFFile,i,0,gffindex);
			WriteEnhancerLine(prs,GFFFile,i,0,prs.refseq[i].AllExperiments_IntBins[0].ClusteredEnrichedBins[j][0],gffindex);
			++gffindex;
		}
	
	}
}

}
void VisualiseInteractions::WriteBEDFiles(PromoterClass& prs,ofstream& BEDFile){
int ColorIndex = 1;


InitialiseColorMatrix();

BEDFile << "track name="  << "Interactions" << "description=" << " BEDformat" << " visibility=2 itemRgb=" << "On" << endl;
for(int i=0;i<prs.NofPromoters;++i){
	for(int j = 0; j < prs.refseq[i].AllExperiments_IntBins[0].ClusteredEnrichedBins.size();++j){
		if( prs.refseq[i].AllExperiments_IntBins[0].ClusteredEnrichedBins[j].size() >= 1){
//			ColorIndex=DetermineColor(prs,i,j,b);
//			if(ColorIndex!=-1){ // If there is an interaction
			WriteBEDLine(prs,BEDFile,i,0,prs.refseq[i].AllExperiments_IntBins[0].ClusteredEnrichedBins[j][0],ColorIndex);
		}
	}
}

}
void VisualiseInteractions::WritePromoterLine(PromoterClass &prs,ofstream& gfile,int prind,int isind,int gffindex){
	gfile << prs.refseq[prind].chr << '\t' << "Interactions" << '\t'<< "Promoter" << '\t';
	gfile << prs.refseq[prind].isoformpromotercoords[isind]-(BinSize/2) << '\t' << prs.refseq[prind].isoformpromotercoords[isind]+(BinSize/2) << '\t';
	gfile << "." << '\t' << "." << '\t' << "." << '\t' << prs.refseq[prind].TranscriptName[isind] << "_" << gffindex << endl;
}
void VisualiseInteractions::WriteEnhancerLine(PromoterClass &prs,ofstream& gfile,int prind,int isind,int binindex,int gffindex){
	int diff;
	gfile << prs.refseq[prind].chr << '\t' << "Interactions" << '\t' << "Enhancer" << '\t';
	if(prs.refseq[prind].strand=="+")
		diff=(binindex-NofInteractorBins)*BinSize;
	else
		diff=(-1*((binindex-NofInteractorBins)*BinSize));

	gfile << ((prs.refseq[prind].isoformpromotercoords[isind]+diff)) << '\t' << ((prs.refseq[prind].isoformpromotercoords[isind]+diff)+BinSize) << '\t';
	gfile << "." << '\t' << "." << '\t' << "." << '\t' << prs.refseq[prind].TranscriptName[isind] << "_" << gffindex <<  endl;
	
}

int VisualiseInteractions::DetermineColor(PromoterClass& prs,int prindex, int isoformindex,int binindex){

// If there are more than one experiment
	if(prs.refseq[prindex].AllExperiments_IntBins[0].Interactor[isoformindex][binindex] == 1 && 
	   prs.refseq[prindex].AllExperiments_IntBins[1].Interactor[isoformindex][binindex] == 1 && 
	   prs.refseq[prindex].AllExperiments_IntBins[2].Interactor[isoformindex][binindex] == 1)
		return -1;

	if(prs.refseq[prindex].AllExperiments_IntBins[0].Interactor[isoformindex][binindex]<SignificanceThreshold && 
	   prs.refseq[prindex].AllExperiments_IntBins[1].Interactor[isoformindex][binindex]<SignificanceThreshold){
		if(prs.refseq[prindex].AllExperiments_IntBins[2].Interactor[isoformindex][binindex]<SignificanceThreshold)
			return 6;
		else
			return 3;
	}
	else{
		if(prs.refseq[prindex].AllExperiments_IntBins[2].Interactor[isoformindex][binindex]<SignificanceThreshold)
			return 4;
		else
			return 0;
	}
	if(prs.refseq[prindex].AllExperiments_IntBins[1].Interactor[isoformindex][binindex]<SignificanceThreshold && 
	   prs.refseq[prindex].AllExperiments_IntBins[2].Interactor[isoformindex][binindex]<SignificanceThreshold)
		return 5;
	else
		return 1;
	if(prs.refseq[prindex].AllExperiments_IntBins[3].Interactor[isoformindex][binindex]<SignificanceThreshold)
		return 2;
	
}




void VisualiseInteractions::WriteBEDLine(PromoterClass& prs,ofstream& bfile,int prind,int isind, int binindex, int colorindex){
/*
chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
chr7    127472363  127473530  Pos2  0  +  127472363  127473530  255,0,0
chr7    127473530  127474697  Pos3  0  +  127473530  127474697  255,0,0
chr7    127474697  127475864  Pos4  0  +  127474697  127475864  255,0,0
chr7    127475864  127477031  Neg1  0  -  127475864  127477031  0,0,255
chr7    127477031  127478198  Neg2  0  -  127477031  127478198  0,0,255
chr7    127478198  127479365  Neg3  0  -  127478198  127479365  0,0,255
chr7    127479365  127480532  Pos5  0  +  127479365  127480532  255,0,0
chr7    127480532  127481699  Neg4  0  -  127480532  127481699  0,0,255
*/
	
	int diff;
	bfile << prs.refseq[prind].chr << '\t';
	if(prs.refseq[prind].strand=="+")
		diff=(binindex-NofInteractorBins)*BinSize;
	else
		diff=(-1*((binindex-NofInteractorBins)*BinSize));

	bfile << ((prs.refseq[prind].isoformpromotercoords[isind]+diff)) << '\t' << ((prs.refseq[prind].isoformpromotercoords[isind]+diff)+BinSize) << '\t';
	bfile << prs.refseq[prind].TranscriptName[isind] << '\t' << 900 << '\t' << "." << '\t';
	bfile << (prs.refseq[prind].isoformpromotercoords[isind]+diff) << '\t' << ((prs.refseq[prind].isoformpromotercoords[isind]+diff)+BinSize) << '\t';
	bfile << Colors[colorindex].color[0] << "," << Colors[colorindex].color[1] << "," << Colors[colorindex].color[2] << endl;

}