#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#ifndef ARGUMENT_HH
#define ARGUMENT_HH

namespace Barcode{
struct Argument
{
	unsigned overlap;
	unsigned kmerll;
	unsigned kmerInterval;
	float lowFreCut;
	unsigned min;
	unsigned max;
//	unsigned minCov;
//	unsigned avgCov;
	unsigned thrd_n;
	float readSup;
	std::string endsFile;
	std::string outPrefix;
//	std::string readsFile1;
//	std::string readsFile2;
	std::string readsList;

	Argument();
	void fromArg(int argc,char *argv[]);
};
}

#endif
