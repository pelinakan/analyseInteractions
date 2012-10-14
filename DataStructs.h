vector < string > AbNames;
struct EnrichmentperTSS{
	vector <double> peakbins;
	vector <int> nofpeakbins;
	vector <double> nonpeakbins;
	vector <int> nofnonpeakbins;
};
struct AllPeakStruct{
	vector < vector < int > > PeakBins; 
};
struct InteractionStruct{
	vector < vector < double > > interactorbins; // the number of reads in the bins
	vector < vector < double > > Interactor; // If a bin has a distal interactor
	vector < vector < double > > ClusteredEnrichedBins;
};

struct FeatureStruct{
		vector <int> ProbeIDs;
		string FeatureType; //Promoter, NegCtrl, etc..
		string chr;
		vector < InteractionStruct > AllExperiments_IntBins;
		vector < AllPeakStruct > AllPeaks_PeakBins;
//Used for Promoters Only
		string RefSeqName; 
		vector <string> TranscriptName;
		vector <int> isoformpromotercoords; 
		int TSS; //Transcription Start Site
		int gene_end; 
		string strand;
		double *expression;
		EnrichmentperTSS EnrPerTSS;

//Used for NegativeControls Only
		string type; //genic or intergenic
		int start; 
		int end; 
		int midpoint;
};
struct PairStruct{
	string chr_1;
	string chr_2;
	int startcoord;
	int endcoord;
};
/*
struct PromoterStruct{
		vector <int> ProbeIDs;
		string RefSeqName;
		string chr;
		vector <int> isoformpromotercoords;
		int TSS; //Transcription Start Site
		int gene_end; //Gene End Coordinate
		string strand;
		int *interactorbins; // the number of reads in the bins
		int *PeakBins; 
		double *expression;
		EnrichmentperTSS EnrPerTSS;
};
*/
/*
struct NegCtrlStruct{
		vector <int> ProbeIDs;
		string chr;
		string type; //genic or intergenic
		int start; //Transcription Start Site
		int end; //Gene End Coordinate
		int midpoint;
		int *interactorbins; // first element is the index of the bin, the second is the frequency
		int *PeakBins;
};
*/
