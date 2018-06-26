#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <gzstream.h>

#include <cstdlib>

using namespace std;

class fileHandle {
	public:
		fileHandle(std::string file) : file_(file) {
			if(file.substr(file.size()-3) == ".gz")
				FTYPE_ = GZS;
			else FTYPE_ = FS;
		};
		~fileHandle() { }

		void fileOpen() {
			if(FTYPE_ == GZS)  gzs.open(file_.c_str());
			else fs.open(file_.c_str());
		}
		void fileClose() {
			if(FTYPE_ == GZS) gzs.close();
			else fs.close();
		}
		std::string readLine () {
			std::string line;
			if(FTYPE_ == GZS) getline(gzs, line);
			else getline(fs, line);

			return line;
		}
		bool ifend() {
			if(FTYPE_ == GZS)  return gzs.eof();
			else return fs.eof();
		}
	private:
		std::string file_;
		enum FTYPE {GZS, FS} FTYPE_;
		igzstream gzs;
		ifstream fs;

};

void usage () {
	std::cout << '\n' << "Usage: compare vcf file" << '\n';
	std::cout << '\n' << "Input: <benchmark vcf[.gz]> <case vcf[.gz]> [bed[.gz]]" << '\n' << std::endl;
	exit(-1);
}


void readFile (map<std::string,std::vector<std::string>>& vcfInfo, const std::string& file, long long& total) {
	fileHandle fH(file);
	fH.fileOpen();
	while(!fH.ifend()) {
		std::string line = fH.readLine();
		if(line.empty()) continue;
		if(line[0] == '#')  continue;
		istringstream record(line);
		std::vector<std::string> data;
		std::string info;
		while(record >> info) {
			data.push_back(info);
		}
		std::string chro = data[0];
		std::string pos = data[1];
		std::string alt = data[4];
		std::string score = data[5];
		std::vector<std::string> tmp{alt, score};

		vcfInfo.insert({chro+"|"+pos+"|"+alt, tmp});
		//vcfInfo.insert({chro+"|"+pos, tmp});

		total++;
	}
	fH.fileClose();
}

void diff(map<std::string, std::vector<std::string>>& bench_, map<std::string, std::vector<std::string>>& case_, long long& consistence) {
	for (auto it = bench_.begin(); it != bench_.end(); it++) {
		std::string key = it->first;
		if (case_.count(key)) {
			consistence++;
		}
	}
}


int main (int argc, char* argv[]) {
	bool bedMode = false;
	std::string benchVcf, caseVcf, bed;
	if(argc == 4) {
		bedMode = true;   bed = argv[3];
	}
	else if(argc != 3)     usage();
	benchVcf = argv[1];    caseVcf = argv[2];
	long long benchTotal=0, caseTotal=0, consistence=0;
	map<std::string, std::vector<std::string>> vcfInfo1, vcfInfo2;
	readFile(vcfInfo1, benchVcf, benchTotal);
	readFile(vcfInfo2, caseVcf, caseTotal);
	diff(vcfInfo1, vcfInfo2, consistence);
	long long FN = benchTotal - consistence;
	long long FP = caseTotal - consistence;
	double precision = (double)consistence/caseTotal;
	double recall = (double)consistence/benchTotal;
	double Fmeasure = (double)2*precision*recall/(precision+recall);
	std::cout << "bench\tcase\tconsistence\tpresion\trecall\tFmasure\n";
	std::cout << benchTotal << '\t' << caseTotal << '\t' \
	          << consistence << '\t' << precision << '\t' \
		  << recall << '\t' << Fmeasure << std::endl;

	return 0;
}

