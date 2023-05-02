void find_bad_files() // run this in container at /gpfs01/star/pwg/youqi/run15/result with
					  // root -l ../analysis_code/find_bad_files.C
{
	string filelistname = "/gpfs01/star/pwg/youqi/run15/result/allfilelist.list";
	const char *c = filelistname.c_str();
	FILE *filelist = fopen(c, "r");
	char file[256];
	while (fgets(file, 256, filelist) != NULL)
	{
		string filestr;
		filestr += file;
		filestr.erase(remove(filestr.begin(), filestr.end(), '\n'), filestr.end());
		TFile f(filestr.c_str());
        if (f.IsZombie())
		{
			cout << filestr << endl;
			continue;
		}
	}
}