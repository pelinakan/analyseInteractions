class RESitesClass{
	friend class PromoterClass;
public :
	int span;
	int window;

	vector <string> chr_names;
	vector <int> chr_offsets; //offset bit to start at the right chr
	vector < int > chrstarts;
	vector < int > chrends;
	int MappTextBinaryFileSize;
	void InitialiseVars(void);
	void CreateIndexFile_RESites(void);
	void WriteRESitesText_toBinaryFile(void);
	int GetRESitesCount(string, int);
private:
	void FindBinaryPos(int, int, int, int&, int&);


};

void RESitesClass::CreateIndexFile_RESites(void){
	ifstream RESitesf("mm9_GATCpos.txt");

	ofstream REindexf("MM9.GATCPos.index.txt");// Format : chr '\t' start '\t' offset
	ofstream REbinaryf("MM9.GATCPos.bin", ios::binary); // start end mappability
	ofstream REchrbounds("MM9.chrLengths.txt");

	string chrname,chrp,temp;
	int pos, chrstart,binst, binend; 
	int span = 1000000; // Window size
	int count = 0; // how many windows per chunk
	int chrcount = 0;
	int offset;
	vector < int > posvector;

	//For indexing
	int icount = 0;
	int chunksize = 1000000; // 1 megabase chuncks

	bool f=0;
	RESitesf >> chrp >> pos >> temp >> temp; // Read the first line
	chrname = chrp; // get the chr name outside the loop
	chrstart = pos; // chromosome start
	
	while(!(RESitesf.eof())){
		int binstart = pos;
		int binend = binstart + span;
		while(chrname == chrp && (pos >= binstart && pos <= binend)){
			posvector.push_back(pos);
			++count;
			RESitesf >> chrp >> pos >> temp >> temp; // Read the first line
			if(RESitesf.eof()){
				f = 1; //End of file
				break;
			}
		}
		int offset = REbinaryf.tellp();
		for(int i=0; i < posvector.size();++i)
			REbinaryf.write((char *)(&(posvector[i])), sizeof((posvector[i])));
		
		REindexf << chrname << '\t' << binstart << '\t' << binend << '\t'
					<< span << '\t' << posvector.size() << '\t' << offset << endl;

		if((chrname != chrp) || f ){			
			REchrbounds << chrname << '\t' << chrstart << '\t' << posvector.back() << endl;
			chrname = chrp;
			chrstart = pos;
			cout << chrname <<  endl;
		}
		posvector.clear();
		count = 0;

	}
}
void RESitesClass::WriteRESitesText_toBinaryFile(void){
string temp, chr, chrtemp;

int start, end, count, offset, fsize;
	ifstream MappF("MM9.GATCPos.index.txt");
	ofstream MappFBin("MM9.GATC.offsets.binary", ios::binary);
// Mapptext_mm9.binary file contains start, end count and offset information

	fsize = MappFBin.tellp();
	chr_offsets.push_back(fsize);
	MappF >> chr >> start >> end >> temp >> count >> offset;
	chr_names.push_back(chr);
	getline(MappF,temp);
	MappFBin.write((char *)(&start), sizeof(start));
	MappFBin.write((char *)(&end), sizeof(end));
	MappFBin.write((char *)(&offset),sizeof(offset));
	MappFBin.write((char *)(&count), sizeof(count));

	do{
		do{
			MappF >> chrtemp >> start >> end >> temp >> count >> offset;
			getline(MappF,temp);
			MappFBin.write((char *)(&start), sizeof(start));
			MappFBin.write((char *)(&end), sizeof(end));
			MappFBin.write((char *)(&offset),sizeof(offset));
			MappFBin.write((char *)(&count), sizeof(count));
		}while(chr == chrtemp && (!MappF.eof()));
		chr_names.push_back(chrtemp);
		fsize = MappFBin.tellp();
		chr_offsets.push_back(fsize);
		chr = chrtemp;
	}while(!MappF.eof());
	MappFBin.close();

	ofstream ofile("GATC.mm9.chr_offsets.txt");
	for(int i=0;i<chr_names.size();++i)
		ofile << chr_names[i] << '\t' << chr_offsets[i] << endl;
}
void RESitesClass::InitialiseVars(void){
	span = 200;
	window = 200000;

	ifstream ifile("/bubo/proj/b2011029/bin/Detect_Interactions/supportingFiles/GATC.mm9.chr_offsets.txt");
	
	ifstream Mapptextbin("/bubo/proj/b2011029/bin/Detect_Interactions/supportingFiles/MM9.GATC.offsets.binary", ios::binary);
	Mapptextbin.seekg(0, ios::end);
	MappTextBinaryFileSize = Mapptextbin.tellg();
	Mapptextbin.close();

	string chrname;
	int chroffset;
	do{
		ifile >> chrname >> chroffset;
		chr_names.push_back(chrname);
		chr_offsets.push_back(chroffset);
	}while(!ifile.eof());
	ifile.close();

	ifstream chrbounds("/bubo/proj/b2011029/bin/Detect_Interactions/supportingFiles/MM9.chrLengths.txt");
	int st, end;
	chrstarts.resize(chr_names.size());
	chrends.resize(chr_names.size());
	do{
		chrbounds >> chrname >> st >> end;
		for (int i = 0; i < chr_names.size();++i){
			if(chrname == chr_names[i]){
				chrstarts[i] = st;
				chrends[i] = end;
				break;
			}
		}
	}while(!chrbounds.eof());
	chrbounds.close();




}
void  RESitesClass::FindBinaryPos(int chroffset, int chroffset_next, int coord, int &fileoffset, int &bitcount){

ifstream file2 ("/bubo/proj/b2011029/bin/Detect_Interactions/supportingFiles/MM9.GATC.offsets.binary",ios::binary);
// GCtext_mm9.binary file contains start, end count and offset information

int start, end, count, offset;

    file2.seekg(chroffset, ios::beg); // Get to the correct chromosome position 
	do{
		file2.read((char *)(&start),sizeof(start));
		file2.read((char *)(&end),sizeof(end));
		file2.read((char *)(&offset),sizeof(offset));
		file2.read((char *)(&count),sizeof(count));
		
		if(( start <= coord ) && ( end >= coord )){
				fileoffset = offset; // this is the offset
				bitcount = count;
				if ((coord + BinSize) > end ){
					file2.read((char *)(&start),sizeof(start));
					file2.read((char *)(&end),sizeof(end));
					file2.read((char *)(&offset),sizeof(offset));
					file2.read((char *)(&count),sizeof(count));
					
					bitcount += count;
				}
				break;
		}
	}while (!(file2.eof()));

	file2.close();

}
int RESitesClass::GetRESitesCount(string chr, int coord){

ifstream file ("/bubo/proj/b2011029/bin/Detect_Interactions/supportingFiles/MM9.GATCPos.bin", ios::in | ios::binary | ios::ate);

ifstream::pos_type fileSize;

int bitcount = 0;
int REcount = 0;
if(file.is_open())
{
	int offset;
	int chroffset,chroffset_next;
	int binstart = coord;
	int binend = coord + BinSize;

	for(int i=0;i<chr_names.size();++i){
		if(chr_names[i] == chr){
			if ((binend <= chrstarts[i] ) || binstart >= chrends[i] )
				return 0;
			if (binstart <= chrstarts[i] && (binend >= chrstarts[i]))
				binstart = chrstarts[i];

			chroffset = chr_offsets[i];
			if(i<chr_names.size())
				chroffset_next = chr_offsets[i+1];
			else
				chroffset_next = MappTextBinaryFileSize;
			break;
		}
	}	
	FindBinaryPos(chroffset, chroffset_next, binstart, offset, bitcount);

	int REpos = 0;

	file.seekg(offset, ios::beg);
	for (int i = 0; i < bitcount ; ++i){
		if (file.eof())
			break;
			file.read((char *)(&REpos),sizeof(REpos));
		while (REpos >= binstart && REpos <= binend ){
			if (file.eof())
				break;
				file.read((char *)(&REpos),sizeof(REpos));
			++REcount;
		}
		if(REcount > 0)
			break;
	}
	file.close();
}
	return REcount;

}
