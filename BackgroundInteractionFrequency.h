/*
Background interaction frequency is calculated as 
mean interaction frequency per bin for all promoters
*/
class DetermineBackgroundLevels{
public:
	double* Distance;
	double* log2Distance;
	double* mean;
	double* stdev;
	void InitialiseVars(void);
	void CalculateMeanandStd(PromoterClass&,NegCtrlClass&,string,int);
	void CalculateMeanandStdRegress(PromoterClass&,string,int);

	void ClearBackgroundLevels();

private:
	void LinearRegression(int,double*,double*,string);

};
void DetermineBackgroundLevels::ClearBackgroundLevels(){
	delete []mean;
	delete []stdev;
}
void DetermineBackgroundLevels::InitialiseVars(void){

	Distance = new double[NofInteractorBins];
	log2Distance = new double[NofInteractorBins];
	int i;
	Distance[0]=1;
	for(i=1;i<NofInteractorBins;++i){
		Distance[i]=BinSize*i;
		log2Distance[i]=(log(double(BinSize*i))/log(2.0));
	}

	mean = new double[NumberofBins];
	stdev = new double[NumberofBins];
}
void DetermineBackgroundLevels::CalculateMeanandStdRegress(PromoterClass& prs, string whichstream,int whichExperiment){
int i,j,k,offset,index,start,count=0;

double *mean,*stdev,mean_square;

mean = new double[NofInteractorBins];
stdev = new double[NofInteractorBins];

if(whichstream=="upstream"){
	start=1;
	offset=-1;
	index=(NofInteractorBins);
}
else{
	start=NofInteractorBins;
	offset=1;
	index=0;
}
for(i=start;;i++){
	mean[index]=0.0;
	stdev[index]=0.0;
	for(j=0;j<prs.NofPromoters;j++){
		for(k=0;k<prs.refseq[i].isoformpromotercoords.size();++k){
			mean[index]+=prs.refseq[j].AllExperiments_IntBins[whichExperiment].interactorbins[k][i];
			stdev[index]+=pow((double(prs.refseq[j].AllExperiments_IntBins[whichExperiment].interactorbins[k][i])),2);
		}
	}
	mean[index]/=(double(prs.NofPromoters));
	mean_square=stdev[index]/(double(prs.NofPromoters));
	stdev[index]=(mean_square-(pow(mean[index],2)));
	stdev[index]/=(sqrt(double(prs.NofPromoters)));
	mean[index]=(log(mean[index])/log(2.0));
	++count;
	index+=offset;
	
	if(count>NofInteractorBins)
		break;
}	
	LinearRegression(NofInteractorBins,Distance,mean,whichstream);
}
void DetermineBackgroundLevels::LinearRegression(int n, double* x, double* y,string whichstream){
double j=0;
int k=0,start,offset;
double intercept,a,b; //Power law y=a*x^b
Maths::Regression::Linear A(n, x, y);
 
    cout << "    Slope = " << A.getSlope() << endl;
	b=A.getSlope();
	cout << "Intercept = " << A.getIntercept() << endl << endl;
	intercept=A.getIntercept();
    a=pow(2,intercept);
	

    cout << "Regression coefficient = " << A.getCoefficient() << endl;
    cout << endl << "Regression line values" << endl << endl;
	if(whichstream=="upstream"){
		start=NofInteractorBins-1;
		offset=-1;
	}
	else{
		start=0;
		offset=1;
	
	}
    for ( int i = start;; i +=offset)  {	
		k++;
		mean[k]=(a*pow(Distance[i],b));
		cout << "x = " << Distance[i] << "  y = " << mean[k] << endl;
		j++;
		if(j>=NofInteractorBins)
			break;
	}
}


/*
void DetermineBackgroundLevels::CalculateMeanandStd(PromoterClass& prs){
int i,j,nofproms=0;

for(i=0;i<NumberofBins;i++){
	mean_norm[i]=0.0;
	stdev_norm[i]=0.0;
	for(j=0;j<prs.NofPromoters;j++){
		if(prs.refseq[i].interactorbins[NofInteractorBins]!=0){
			mean_norm[i]+=prs.refseq[j].interactorbins[i];
			stdev_norm[i]+=pow((prs.refseq[j].interactorbins[i]),2);
			++nofproms;
		}
	}
	mean_norm[i]/=(double(nofproms));
	stdev_norm[i]/=(double(nofproms));
	stdev_norm[i]=(stdev_norm[i]-(pow(mean_norm[i],2)));
	stdev_norm[i]=sqrt(stdev_norm[i]);
}	
}
*/
void DetermineBackgroundLevels::CalculateMeanandStd(PromoterClass& prs, NegCtrlClass& nc,string BaseFileName,int whichExperiment){
int i,j,k;
vector < int > nofncs;
nofncs.resize(NumberofBins);
string filen;
for(i=0;i<NumberofBins;i++){
	mean[i]=0.0;
	stdev[i]=0.0;
}

for(i=1;i<NumberofBins;i++){
	nofncs[i]=0;
	for(j=0;j<prs.NofPromoters;j++){
		for(k=0;k<prs.refseq[j].isoformpromotercoords.size();k++){
			if(prs.refseq[j].AllExperiments_IntBins[whichExperiment].interactorbins[k][i]<ExcludeOutliers){
				mean[i]+=prs.refseq[j].AllExperiments_IntBins[whichExperiment].interactorbins[k][i];
				stdev[i]+=pow(double((prs.refseq[j].AllExperiments_IntBins[whichExperiment].interactorbins[k][i])),2);
				++nofncs[i];
			}
		}
	}
	mean[i]/=(double(nofncs[i]));
	stdev[i]/=(double(nofncs[i]));
	stdev[i]=(stdev[i]-(pow(mean[i],2)));
	stdev[i]=sqrt(stdev[i]);
}
/*
for(i=1;i<NumberofBins;i++){
	for(j=0;j<nc.NofNegCtrls;j++){
//		if(nc.negctrls[j].interactorbins[0][NofInteractorBins]>1){
			mean[i]+=nc.negctrls[j].interactorbins[0][i];
			stdev[i]+=pow(double((nc.negctrls[j].interactorbins[0][i])),2);
			++nofncs[i];
//		}
	}
	mean[i]/=(double(nofncs[i]));
	stdev[i]/=(double(nofncs[i]));
	stdev[i]=(stdev[i]-(pow(mean[i],2)));
	stdev[i]=sqrt(stdev[i]);
}	
*/
#ifdef UNIX
	filen.append(dirname.c_str());
#endif
filen.append(BaseFileName);
filen.append("BackgroundLevels.txt");

ofstream outf(filen.c_str());
for(i=1;i<NumberofBins;i++)
	outf <<	mean[i] << '\t' <<	stdev[i] << endl;

outf.close();

}