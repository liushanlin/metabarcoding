#include "Argument.hh"

using namespace std;
using namespace Barcode;
Argument::Argument():
	overlap(121),
	kmerInterval(10),
	kmerll(55),
	lowFreCut(0.1),
	max(700),
	min(450),
	readSup(0.1),
	thrd_n(8),
	endsFile(""),
	outPrefix("out"),
	readsList(""){}

void Argument::fromArg(int argc,char *argv[]){
	string Usage="\
Usage: SOAPbarcode [options]\n\
\n\
Options(default):\n\
	-e: ends file\n\
	-r: reads list\n\
	-k: length of kmer(121),permit 127mer\n\
	-l: kmer lower limit(55)\n\
	-v: kmer interval(10)\n\
	-c: cutoff of low frequence kmer(0.1)\n\
	-s: reads support number while solve branches(0.1)\n\
	-x: maximal length of sequence assembly(700)\n\
	-n: minimal length of sequence assembly(450)\n\
	-t: CPU number(8)\n\
	-o: prefix of out put file(out)\n\
	-h: help\n\
	\n\
Report bug to: lujianliang@genomics.org.cn shanlin.liu@ioz.ac.cn\n\
";         //kmerul and kmerInterval don't illustrate yet
	if (argc<2)
	{
		cout<<Usage;
		exit(1);
	}
	int c;
	while ((c=getopt(argc,argv,"e:r:k:l:v:c:x:n:t:o:s:h"))!=-1){
		switch (c)
		{
		case 'e':
			endsFile=optarg;
			break;
		case 'r':
			readsList=optarg;
			break;
		case 'k':
			overlap=atoi(optarg);
			break;
		case 'l':
			kmerll=atoi(optarg);
			break;
		case 'v':
			kmerInterval=atoi(optarg);
			break;
		case 'c':
			lowFreCut=atof(optarg);
			break;
		case 'x':
			max=atoi(optarg);
			break;
		case 'n':
			min=atoi(optarg);
			break;
		case 't':
			thrd_n=atoi(optarg);
			break;
		case 'o':
			outPrefix=optarg;
			break;
		case 's':
			readSup=atof(optarg);
			break;
		case 'h':
			cout<<Usage;
			exit(1);
		case '?':
			cout<<Usage;
			exit(1);
		}
	}
	if (overlap%2==0 || overlap>127 || overlap<25)
	{
		cout<<Usage;
		exit(1);
	}
	if (kmerll<25){
		cout<<Usage;
		exit(1);
	}
	if (kmerInterval%2==1){
		cout<<Usage;
		exit(1);
	}
	if (endsFile=="" || readsList=="")
	{
		cout<<Usage;
		exit(1);
	}
	if (readSup<0)
	{
		cout<<Usage;
		exit(1);
	}
/*	if (minCov<1 || minCov>avgCov)
	{
		cout<<Usage;
		exit(1);
	}
	if (avgCov<minCov || avgCov>63)
	{
		cout<<Usage;
		exit(1);
	}*/
	if (thrd_n>300)
	{
		cout<<Usage;
		exit(1);
	}
	if (lowFreCut<0 || lowFreCut>1)
	{
		cout<<Usage;
		exit(1);
	}
}

