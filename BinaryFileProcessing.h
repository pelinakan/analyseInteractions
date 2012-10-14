class MappabilityClass{

public :
public :
	int span;
	int window;
	int NumberofBits;
	double lowerLimit;
	double dataRange;

	vector <string> chr_names;
	vector <int> chr_offsets; //offset bit to start at the right chr
	int MappTextBinaryFileSize;
	void InitialiseVars(void);
	void CreateIndexFile_Mappability(void);
	void WriteMappabilityText_toBinaryFile(void);
	double GetMappability(string, int);
private:
	void FindBinaryPos(int, int, int, int&, int&);
	int CalculateNumberofBitsToRead(ifstream&, int&, int, int&);


};

void MappabilityClass::CreateIndexFile_Mappability(void){
	ifstream mappf("crgMapabilityAlign36mer.bedgraph");

	ofstream mappindexf("MM9.Mappability.36mer.index.txt");// Format : chr '\t' start '\t' offset
	ofstream mappbinaryf("MM9.Mappability.36mer.bin", ios::binary); // start end mappability

	string chrname,chrp, chrn;
	int startp,startn, endp, endn, chrstart,binst, binend; 
	double mp,mn,nextbinm;
	int span = 200; // Window size
	int count = 500; // how many windows per chunk
	int offset;
	vector < int > starts;
	vector <int > ends;
	vector < double > mapp;
	vector < double > tm;

	//For indexing
	int icount = 0;
	int chunksize = 1000000; // 1 megabase chuncks

	bool f=0;
	mappf >> chrp >> startp >> endp >> mp; // Read the first line
	while(!(mappf.eof())){
		chrname = chrp; // get the chr name outside the loop
		chrstart = startp; // chromosome start
		int binstart = chrstart;
		int binend = binstart + span;
		cout << chrname << "  " << chrstart << "  " << endl;
		while(chrname == chrp &&(startp >= binstart && startp <=binend)){
			if (endp <= binend) // it belongs to the next bin also
				tm.push_back((mp*(endp-startp)));				
			else{
				tm.push_back(mp*(binend-startp));
				//check how many bins it contains
				int remainder = endp - binend;
				int rbins = remainder/span;
				double meanm = 0.0;
				if (rbins == 0){
					for(int it=0; it<tm.size();++it)
						meanm +=tm[it];
					mapp.push_back(meanm/(double(span)));
					starts.push_back(binstart);
					ends.push_back(binend);
					binstart +=span;
					binend +=span;
				}
				else{ // if it contains more bins
					for(int it=0; it<tm.size();++it)
						meanm +=tm[it];
					mapp.push_back(meanm/(double(span)));
					starts.push_back(binstart);
					ends.push_back(binend);
					binstart +=span;
					binend +=span;
					for (int b = 1; b <= rbins; ++b){
						mapp.push_back(mp);
						starts.push_back(binstart);
						ends.push_back(binend);
						binstart +=span;
						binend +=span;
					}
				}
				tm.clear();
				tm.push_back(mp*(endp-binstart));
			}
			int prevend = endp;
			mappf >> chrp >> startp >> endp >> mp; // Read the first line
			if(mappf.eof())
				break;
			if(chrname == chrp && (startp != prevend)){
				int gapsize = startp - prevend;
				if( binend >= startp) // If the next line is within the current bin
					tm.push_back(0);
				else{
					//First finish the current bin
					tm.push_back(0);
					double meanm = 0.0;
					for(int it=0; it<tm.size();++it)
						meanm +=tm[it];
					mapp.push_back(meanm/(double(span)));
					starts.push_back(binstart);
					ends.push_back(binend);
					int gapbins = (startp - binend)/span;
					binstart+=span;
					binend+=span;
					// If there are more bins in the gap
					for (int g = 1; g <= gapbins; ++g){
						mapp.push_back(0);		
						starts.push_back(binstart);
						ends.push_back(binend);
						binstart+=span;
						binend+=span;
					}
					tm.clear();
					if(startp >= binstart)
						tm.push_back(0);
				}
			}
			if(mapp.size() >= count){
				int offset = mappbinaryf.tellp();
				for(int i = 0; i < mapp.size(); ++i){
//					mappbinaryf.write((char *)(&(mapp[i])), sizeof((mapp[i])));
//					int offset2 = mappbinaryf.tellp();
					int compressed = 127*(mapp[i]);
//					mappbinaryf.write((char *)(&(compressed)), sizeof((compressed)));
//					int offset3 = mappbinaryf.tellp();
					char ascii =char(compressed);
					mappbinaryf.write((char *)(&(ascii)), sizeof((ascii)));
//					int offset4 = mappbinaryf.tellp();
//					cout << offset4 << endl;
//					double value_decoded = lowerLimit+(dataRange*((double(ascii))/127.0));
//					cout << mapp[i] << "  " << compressed << "   " << value_decoded << endl;
				}
				mappindexf << chrname << '\t' << starts[0] << '\t' << ends[(mapp.size()-1)] << '\t'
						   << span << '\t' << mapp.size() << '\t' << offset << endl;
				starts.clear();
				ends.clear();
				mapp.clear();
			}
		}
	}
}



void MappabilityClass::WriteMappabilityText_toBinaryFile(void){
string temp, chr, chrtemp;

int start, end, count, offset, fsize;
	ifstream MappF("MM9.Mappability.36mer.index.txt");
	ofstream MappFBin("Mappability.36mers.offsets_mm9.binary", ios::binary);
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

	ofstream ofile("Mappability.36mers.mm9.chr_offsets.txt");
	for(int i=0;i<chr_names.size();++i)
		ofile << chr_names[i] << '\t' << chr_offsets[i] << endl;
}
void MappabilityClass::InitialiseVars(void){
	span = 200;
	window = 200000;
	NumberofBits = window/span;
	lowerLimit = 0.0;
	dataRange = 1.0;

	ifstream ifile("Mappability.36mers.mm9.chr_offsets.txt");
	
	ifstream Mapptextbin("Mappability.36mers.offsets_mm9.binary", ios::binary);
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

}
void  MappabilityClass::FindBinaryPos(int chroffset, int chroffset_next, int coord, int &fileoffset, int &bitcount){

ifstream file2 ("Mappability.36mers.offsets_mm9.binary",ios::binary);
// GCtext_mm9.binary file contains start, end count and offset information

int start, end, count, offset, NumberofBitsRead;

    file2.seekg(chroffset, ios::beg); // Get to the correct chromosome position

	do{
		file2.read((char *)(&start),sizeof(start));
		file2.read((char *)(&end),sizeof(end));
		file2.read((char *)(&offset),sizeof(offset));
		file2.read((char *)(&count),sizeof(count));
		
		if(( start <= coord ) && ( end >= coord )){
			if(offset != 0)
				fileoffset = (offset + ((coord - start)/span)); // this is the offset
			else
				fileoffset = offset;
             //Check if there is any gaps, determine the bitcount
			 NumberofBitsRead = ((coord - start)/span);
			 bitcount = CalculateNumberofBitsToRead(file2, NumberofBitsRead, coord, end);
			 break;
		}
	}while (!(file2.eof()));

	file2.close();

}
double MappabilityClass::GetMappability(string chr, int coord){

ifstream file ("MM9.Mappability.36mer.bin", ios::in | ios::binary | ios::ate);

ifstream::pos_type fileSize;

char* data;
double Mappmean=0;
int bitcount = 0;
if(file.is_open())
{
	int offset;
	int chroffset,chroffset_next;

	for(int i=0;i<chr_names.size();++i){
		if(chr_names[i] == chr){
			chroffset = chr_offsets[i];
			if(i<chr_names.size())
				chroffset_next = chr_offsets[i+1];
			else
				chroffset_next = MappTextBinaryFileSize;
			break;
		}
	}
	
	FindBinaryPos(chroffset, chroffset_next, coord, offset, bitcount);

	char* data;
	if(bitcount > 0){
		data = new char[bitcount];
		file.seekg(offset, ios::beg);

		if(!file.read(data, bitcount))
            cout << "fail to read" << endl;
		for (int i = 0; i < bitcount ; ++i){
			if ( data[i] <128 ){
				double value = lowerLimit+(dataRange*((double(data[i]))/127.0));
				Mappmean+=value;
				cout << chr << '\t' << (coord + (i*span) + 1 ) << '\t' << value << endl;
			}
		}
		Mappmean/=bitcount;
		file.close();
	}
/*
	if(bitcount > 0){
		file.seekg(offset, ios::beg);

		for (int i = 0; i < bitcount ; ++i){
			double value;
			file.read((char *)(&value),sizeof(value));
			cout << chr << '\t' << (coord + (i*span) + 1 ) << '\t' << value << endl;
			Mappmean+=value;
		}
		Mappmean/=bitcount;
		file.close();
	}
*/
}
	return Mappmean;

}

int MappabilityClass::CalculateNumberofBitsToRead(ifstream &file, int &bitcount, int coord, int &end_prev){

int start, start_next, end, count, offset;

file.read((char *)(&start),sizeof(start)); 
file.read((char *)(&end),sizeof(end));
file.read((char *)(&offset),sizeof(offset));
file.read((char *)(&count),sizeof(count));

if ( end_prev !=  start ) // there is a gap, you cannot read more
	return bitcount; // Return to the number of bits that can be read!
else{
	if ((NumberofBits - bitcount) < count ) // if the remaining bits are all contained in the next step
		return NumberofBits;
	else{
		end_prev = end; // You will read more lines ...
		bitcount += count;
		CalculateNumberofBitsToRead(file, bitcount, coord, end_prev);
	}
}
}

