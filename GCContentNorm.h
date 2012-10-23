class GCContent{
public :
	int span;
	int window;
	int NumberofBits;
	double lowerLimit;
	double dataRange;

	vector <string> chr_names;
	vector <int> chr_offsets; //offset bit to start at the right chr
	int GCTextBinaryFileSize;

	void InitialiseVars(void);
	double GetGCContent(string, int);
	void GenerateGC_bedgraph(void);
	void WriteGCText_toBinaryFile(void);
private:
	void FindBinaryPos(int, int, int, int&, int&);
	int CalculateNumberofBitsToRead(ifstream&, int&, int, int&);
};
void GCContent::InitialiseVars(void){
	span = 5; // Resolution of the file
	window = 5000; // GC content of 100 bp 
	NumberofBits = window/span;
	lowerLimit = 0.0;
	dataRange = 100.0;

	ifstream ifile("/bubo/home/h20/pelin/3Cproj/bin/Detect_Interactions/supportingFiles/GCText_mm9_chr_offsets.txt");
	
	ifstream GCtextbin("/bubo/home/h20/pelin/3Cproj/bin/Detect_Interactions/supportingFiles/GCtext_mm9.binary", ios::binary);
	GCtextbin.seekg(0, ios::end);
	GCTextBinaryFileSize = GCtextbin.tellg();
	GCtextbin.close();

	string chrname;
	int chroffset;
	do{
		ifile >> chrname >> chroffset;
		chr_names.push_back(chrname);
		chr_offsets.push_back(chroffset);
	}while(!ifile.eof());

}
void  GCContent::FindBinaryPos(int chroffset, int chroffset_next, int coord, int &fileoffset, int &bitcount){

ifstream file2 ("/bubo/home/h20/pelin/3Cproj/bin/Detect_Interactions/supportingFiles/GCtext_mm9.binary",ios::binary);
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
double GCContent::GetGCContent(string chr, int coord){

ifstream file ("/bubo/home/h20/pelin/3Cproj/bin/Detect_Interactions/supportingFiles/gc5Base.wib", ios::in | ios::binary | ios::ate);

ifstream::pos_type fileSize;

char* data;
double GCmean=0;
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
				chroffset_next = GCTextBinaryFileSize;
			break;
		}
	}
	
	FindBinaryPos(chroffset, chroffset_next, coord, offset, bitcount);

	if(bitcount > 0){
		data = new char[bitcount];
		file.seekg(offset, ios::beg);

		if(!file.read(data, bitcount))
            cout << "fail to read" << endl;
		for (int i = 0; i < bitcount ; ++i){
			if ( data[i] < 128 ){
				double value = lowerLimit+(dataRange*((double)data[i]/127.0));
				GCmean+=value;
//				cout << chr << '\t' << (coord + (i*span) + 1 ) << '\t' << value << endl;
			}
		}
		GCmean/=bitcount;
		file.close();
	}
}
	return GCmean;

}

void GCContent::GenerateGC_bedgraph(void){

	ifstream file ("GCtext_mm9.binary",ios::binary);
	ifstream::pos_type position;

	int i, start, end, offset, count, num_windows,remainder,coordst, coordend;
	int bedgraphwindowsize = 4; // corresponding to 20 bases

	ifstream wibfile ("gc5Base.wib", ios::in | ios::binary | ios::ate);
	char* data;
	
	ofstream bedgraph_file("mm9_GC_20bases.bedgraph");

	for (i = 0; i<(chr_offsets.size()-1); ++i){
//		file.seekg(chr_offsets[i], ios::beg); // Get to the correct chromosome position
		do{
			
			file.read((char *)(&start),sizeof(start));
			file.read((char *)(&end),sizeof(end));
			file.read((char *)(&offset),sizeof(offset));
			file.read((char *)(&count),sizeof(count));
			
			wibfile.seekg(offset, ios::beg); // Get to the right position
			num_windows = count / bedgraphwindowsize;	
			remainder = count % bedgraphwindowsize;			
			coordst = start;
			for (int j = 0; j < num_windows; ++j){
				coordend = coordst + bedgraphwindowsize*span;
				data = new char[bedgraphwindowsize];
				double GCmean = 0.0;
				double value = 0.0;
				if(!wibfile.read(data, bedgraphwindowsize))
					cout << "fail to read" << endl;
				for (int k = 0; k < bedgraphwindowsize ; ++k){
					if ( data[k] < 128 ){
							value = lowerLimit+(dataRange*((double)data[k]/127.0));
							GCmean+=value;
					}
				}
				GCmean/=bedgraphwindowsize;
				bedgraph_file << chr_names[i] << '\t' << ( coordst + 1 ) << '\t' << ( coordend + 1 ) << '\t' << GCmean << endl;
				coordst += bedgraphwindowsize*span;
			}
			if(remainder > 0 ){
				data = new char[remainder];
				double GCmean = 0.0;
				double value = 0.0;
				if(!wibfile.read(data, remainder))
					cout << "fail to read" << endl;
				for (int k = 0; k < remainder ; ++k){
					if ( data[k] < 128 ){
							value = lowerLimit+(dataRange*((double)data[k]/127.0));
							GCmean+=value;
					}
				}
				GCmean/=remainder;
				bedgraph_file << chr_names[i] << '\t' << ( coordst + 1 ) << '\t' << ( coordst + 1 + (span*remainder) ) << '\t' << GCmean << endl;
			}

			position = file.tellg();
		}while(position < chr_offsets[i+1]);
		cout << chr_names[i] << "   finished" << endl;
	}
	// Do it for the last chromosome
	file.seekg(chr_offsets[i], ios::beg); // Get to the correct chromosome position
	do{
		file.read((char *)(&start),sizeof(start));
		file.read((char *)(&end),sizeof(end));
		file.read((char *)(&offset),sizeof(offset));
		file.read((char *)(&count),sizeof(count));
		file.seekg(offset, ios::beg); // Get to the right position
		num_windows = count / bedgraphwindowsize; 
		for (int j = 0; j < num_windows; ++j){
			data = new char[bedgraphwindowsize];
			double GCmean = 0.0;
			double value = 0.0;
			if(!file.read(data, bedgraphwindowsize))
				cout << "fail to read" << endl;
			for (int k = 0; k < bedgraphwindowsize ; ++k){
				if ( data[k] < 128 ){
					double value = lowerLimit+(dataRange*((double)data[k]/127.0));
					GCmean+=value;
				}
			}
			GCmean/=bedgraphwindowsize;
			cout << chr_names[i] << '\t' << (start + (j*span*bedgraphwindowsize) + 1 ) << '\t' << (start + (j*span) + 1 + (bedgraphwindowsize*span)) << '\t' << value << endl;
		}
		position = file.tellg();
	}while(offset < GCTextBinaryFileSize);

}

int GCContent::CalculateNumberofBitsToRead(ifstream &file, int &bitcount, int coord, int &end_prev){

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
void GCContent::WriteGCText_toBinaryFile(void){
string temp, chr, chrtemp;

int start, end, count, offset, fsize;
	ifstream GC5Base("gc5Base_mm9.txt");
	ofstream GCtextbin("GCtext_mm9.binary", ios::binary);
// GCtext_mm9.binary file contains start, end count and offset information

	fsize = GCtextbin.tellp();
	chr_offsets.push_back(fsize);
	GC5Base >> temp >> chr >> start >> end >> temp >> temp >> count >> offset;
	chr_names.push_back(chr);
	getline(GC5Base,temp);
	GCtextbin.write((char *)(&start), sizeof(start));
	GCtextbin.write((char *)(&end), sizeof(end));
	GCtextbin.write((char *)(&offset),sizeof(offset));
	GCtextbin.write((char *)(&count), sizeof(count));

	do{
		do{
			GC5Base >> temp >> chrtemp >> start >> end >> temp >> temp >> count >> offset;
			getline(GC5Base,temp);
			GCtextbin.write((char *)(&start), sizeof(start));
			GCtextbin.write((char *)(&end), sizeof(end));
			GCtextbin.write((char *)(&offset),sizeof(offset));
			GCtextbin.write((char *)(&count), sizeof(count));
		}while(chr == chrtemp && (!GC5Base.eof()));
		chr_names.push_back(chrtemp);
		fsize = GCtextbin.tellp();
		chr_offsets.push_back(fsize);
		chr = chrtemp;
	}while(!GC5Base.eof());
	GCtextbin.close();

	ofstream ofile("GCText_mm9_chr_offsets.txt");
	for(int i=0;i<chr_names.size();++i)
		ofile << chr_names[i] << '\t' << chr_offsets[i] << endl;

/*
	GC5Base >> temp >> chr >> start >> temp >> temp >> temp >> temp >> offset;
	chr_names.push_back(chr);
	chr_startbits.push_back(offset);
	chr_startcoords.push_back(start);
	getline(GC5Base,temp);
	do{
		do{
			GC5Base >> temp >> chrtemp >> start;
			getline(GC5Base,temp);
		}while(chr == chrtemp && (!GC5Base.eof()));

		chr_names.push_back(chrtemp);
		chr_startbits.push_back(offset);
		chr_startcoords.push_back(start);
		chr = chrtemp;
		if(chr!= "chr1")
			break;
	}while(!GC5Base.eof());

*/
}
