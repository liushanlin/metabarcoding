#include "io.hh"
#include<iostream>

namespace Barcode{
	std::istream& openIn(const std::string& filename,std::ifstream& ifs){
		if (filename=="")
			throw std::runtime_error("Input file name is empty!");
//		std::ifstream ifs; 
		ifs.open(filename.c_str());
		if (!ifs)
			throw std::runtime_error("can't open input file: " + filename);
		return ifs;
	}

	std::ostream& openOut(const std::string& filename,std::ofstream& ofs){
		if (filename=="")
			throw std::runtime_error("output file name is empty!");
//		std::ofstream ofs;
		ofs.open(filename.c_str());
		if (!ofs)
			throw std::runtime_error("can't open output file: " + filename);
		return ofs;
	}
	
	std::ostream& openOut2(const std::string& filename,std::ofstream& ofs){
		if (filename=="")
			throw std::runtime_error("output file name is empty!");
//		std::ofstream ofs;
		ofs.open(filename.c_str(),std::ios::app);
		if (!ofs)
			throw std::runtime_error("can't open output file: " + filename);
		return ofs;
	}
}
