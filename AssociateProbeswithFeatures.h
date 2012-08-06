struct NimblegenProbes{
	string chr;
	int start;
	int end;
	int center;
	int annotated;// 0= not annotated, 1=annotated with a promoter, 2= annotated with a TF site, 3= annotated with a negative control
	bool conflicting_annotations; // if a probe is associated with more than one feature
	int *interactorbins; // the number of reads in the bins, initialised after reading the probes

};

class ProbeSet{
public:
	vector <int> ChrRowStartIndexes;
	vector <int> ChrRowEndIndexes;
	vector <string> ChrNames;
	vector <NimblegenProbes> NotAnnProbeClusters; // Probes that are not negative controls but are not annotated either
	vector <NimblegenProbes> MM9Probes;
	void InitialiseData();
	void ReadProbeCoordinates();
	void ClusterandPrintProbes_NotAnnotated();
	bool AssociateReadwithProbes(string,int,int);
	void PrinttheBins(string);
private:
	void GetProbeFeats(stringstream&,NimblegenProbes&);
};


void ProbeSet::InitialiseData(){
	unsigned int i,j;

	for(i=0;i<MM9Probes.size();++i){
		MM9Probes[i].interactorbins = new int[NumberofBins];
		for(j=0;j<NumberofBins;++j)
			MM9Probes[i].interactorbins[j]=0;
	}

}
void ProbeSet::ReadProbeCoordinates(){

#ifdef UNIX
	string filename;
	filename.append("/bubo/proj/b2011029/bin/3CAnalysis/");
	filename.append("Nimblegen_Capture_Probes_MM9.txt");
	ifstream probefile(filename.c_str());
#endif
#ifdef WINDOWS
	ifstream probefile("Nimblegen_Capture_Probes_MM9.txt");
#endif
	long int index=0;
	string pline,chr1,chr2;

	NimblegenProbes tempprobe;

	getline(probefile,pline);
	stringstream probeline ( pline ); 
	GetProbeFeats(probeline,tempprobe);
	chr1=tempprobe.chr;
	ChrRowStartIndexes.push_back(index);
	ChrNames.push_back(chr1);
	tempprobe.annotated=0;
	tempprobe.conflicting_annotations=0;
	MM9Probes.push_back(tempprobe);
	++index;
	do{
		getline(probefile,pline);
		stringstream probeline ( pline ); 
		GetProbeFeats(probeline,tempprobe);
		chr2=tempprobe.chr;
		tempprobe.annotated=0;
		tempprobe.conflicting_annotations=0;
		MM9Probes.push_back(tempprobe);
		++index;
		while(chr1==chr2 && pline!=""){
			getline(probefile,pline);
			stringstream probeline ( pline ); 
			GetProbeFeats(probeline,tempprobe);
			chr2=tempprobe.chr;
			tempprobe.annotated=0;
			tempprobe.conflicting_annotations=0;
			MM9Probes.push_back(tempprobe);
			++index;
		}
		ChrRowEndIndexes.push_back(index-2);
		if(pline!=""){
			ChrRowStartIndexes.push_back(index-1);
			ChrNames.push_back(chr2);
		}
		chr1=tempprobe.chr;
		cout << "Chr Index   " << index+1 << "   " << chr2 << endl;
	}while(!probefile.eof()&& pline!="");
	MM9Probes.pop_back();

	cout << "Probe Coordinates Read " << endl;
}


void ProbeSet::GetProbeFeats(stringstream& line, NimblegenProbes& t){

	string field;
	getline(line,t.chr,'\t');
	getline(line,field,'\t');
	t.start=atoi(field.c_str());
	getline(line,field,'\t');
	t.end=atoi(field.c_str());
	t.center=t.start+((t.end-t.start)/2);
}



void ProbeSet::ClusterandPrintProbes_NotAnnotated(){


	ofstream outf("Probes_notAssociatedwithAnyFeatures.txt");
	unsigned int c,i,j;
	long int diff;
	NimblegenProbes tempprobecluster;

	for(c=0;c<ChrRowStartIndexes.size();++c){
		for(i=ChrRowStartIndexes[c];i<ChrRowEndIndexes[c];++i){
			if(MM9Probes[i].annotated==0){
				j=i+1;
				do{
					diff=MM9Probes[i].center-MM9Probes[j].center;
					if(MM9Probes[j].annotated!=0)
						break;
					++j;
				}while(diff<=BinSize && j<=i);
				tempprobecluster.start=MM9Probes[i].start;
				tempprobecluster.end=MM9Probes[j-1].end;
				tempprobecluster.center=tempprobecluster.start+((tempprobecluster.end-tempprobecluster.start)/2);
				tempprobecluster.chr=ChrNames[c];
				NotAnnProbeClusters.push_back(tempprobecluster);
			}
		}
	}

	for(i=0;i<NotAnnProbeClusters.size();++i)
		outf << NotAnnProbeClusters[i].chr << '\t' <<  NotAnnProbeClusters[i].start << '\t' <<  NotAnnProbeClusters[i].end << endl;

}
bool ProbeSet::AssociateReadwithProbes(string pchr, int readstart,int readend){

int i=0,j=0,k=0,x=0;
long int startsearch,endsearch;
bool onprobe=0;

startsearch=0;endsearch=-1;
for(j=0;j<ChrNames.size();++j){
	if(pchr==ChrNames[j]){
		startsearch=ChrRowStartIndexes[j];
		endsearch=ChrRowEndIndexes[j];
		break;
	}
}
x=0;
for(j=startsearch;j<=endsearch;++j){
	if(((MM9Probes[j].start-padding)<readstart) && ((MM9Probes[j].end+padding)>readstart)){
		onprobe=1;
		return onprobe;
	}
}
return onprobe;
}


void ProbeSet::PrinttheBins(string basefn){

	int i,j;
	string fn;
	fn.append(basefn);
	fn.append("Probe_InteractorBins.txt");

	ofstream outf(fn.c_str());

	for(i=0;i<MM9Probes.size();++i){
		outf << MM9Probes[i].chr << '\t' << MM9Probes[i].center << '\t';
		for(j=0;j<NumberofBins;++j)
			outf << MM9Probes[i].interactorbins[j] << '\t';
		outf << endl;
	}


}