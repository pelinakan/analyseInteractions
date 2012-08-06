void GenerateFileNames_noEnr_Promoters(string& fn1,string& fn2, string& fn3,string& fn4,string ext,string BaseFileName){
	fn1.append(BaseFileName);
	fn1.append(ext);
	fn1.append("_exp");
	fn1.append(".txt");

	fn2.append(BaseFileName);
	fn2.append(ext);
	fn2.append("_notexp");
	fn2.append(".txt");

	fn3.append(BaseFileName);
	fn3.append(ext);
	fn3.append("_exp");
	fn3.append("_nodir.txt");

	fn4.append(BaseFileName);
	fn4.append(ext);
	fn4.append("_notexp");
	fn4.append("_nodir.txt");

}
void GenerateFileNames_noEnr_NegCtrls(string& fn1,string& fn2,string ext,string BaseFileName){
	fn1.append(BaseFileName);
	fn1.append(ext);
	fn1.append(".txt");

	fn2.append(BaseFileName);
	fn2.append(ext);
	fn2.append("_nodir.txt");
}

void GenerateFileNames_Promoters(string& fn1,string& fn2, string& fn3,string& fn4,string& fn5,string abname, bool enrichment,string BaseFileName){

	fn1.append(BaseFileName);
	fn1.append("_proms");
	fn1.append("_exp");
	fn1.append(abname);
	fn1.append(".txt");

	fn2.append(BaseFileName);
	fn2.append("_proms");
	fn2.append("_notexp");
	fn2.append(abname);
	fn2.append(".txt");

	fn3.append(BaseFileName);
	fn3.append("_proms");
	fn3.append("_exp");
	fn3.append("_nodir");
	fn3.append(abname);
	fn3.append(".txt");

	fn4.append(BaseFileName);
	fn4.append("_proms");
	fn4.append("_notexp");
	fn4.append("_nodir");
	fn4.append(abname);
	fn4.append(".txt");

	if(enrichment){
		fn5.append(BaseFileName);
		fn5.append("_proms_");
		fn5.append("Enrichment");
		fn5.append(abname);
		fn5.append(".txt");
	}

}

void GenerateFileNames(string& fn1,string& fn2,string ext,string BaseFileName){
	fn1.append(BaseFileName);
	fn1.append("_negctrls");
	fn1.append(ext);
	fn1.append(".txt");

	fn2.append(BaseFileName);
	fn2.append(ext);
	fn2.append("_nodir.txt");
}	