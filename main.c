#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <stdio.h>
#include <sstream>
#include <boost/unordered_map.hpp>
#include <pthread.h>
#include "Argument.hh"
#include "Kmer.hh"
#include "Hash.hh"
#include "io.hh"

using namespace std;
using namespace Barcode;

typedef unsigned char uchar;
typedef unsigned int indexT;
typedef boost::unordered_map<unsigned,unsigned> unorderMap;
typedef boost::unordered_multimap<unsigned,unsigned> multiMap;
typedef boost::unordered_map<unsigned,uchar> unordermap;
const unsigned MaxRouteNum=200;
const unsigned MaxRouteLen=650;
const unsigned multiMatch=200;
const int UPlimit=7000;

Argument args;
Hash kset;
Kmer WORDFILTER;
vector<ubyte4> readTkmer;

vector<string> end1;
vector<string> end2;

int endNum=0;
int thrd_num=8;
vector<bool> flagBuf(endNum);
vector<uchar> flagResult(endNum);
vector<uchar> thrdSignal(thrd_num+1);
vector<vector<string> > CONTIG;
vector<vector<int> > traceCount;
static pthread_mutex_t mutex;

int tip_c;

//ofstream os1;
//os1.open("hash.lis");

/*
Kmer kmerBuffer[250];
uchar prevBuffer[250];
uchar nextBuffer[250];
*/
/*static Kmer KmerAnd(Kmer kmer1, Kmer kmer2)
{
	kmer1.high &= kmer2.high;
	kmer1.low &= kmer2.low;
	return kmer1;
}*/

void *ckalloc(unsigned long long amount)
{
	void *p;

	if ((p = (void *) calloc( 1, (unsigned long long) amount)) == NULL && amount != 0)	{
		printf("Ran out of memory while applying %lldbytes\n",amount);
                printf("There may be errors as follows:\n");
                printf("1) Not enough memory.\n");
                printf("2) The ARRAY may be overrode.\n");
                printf("3) The wild pointers.\n");
		fflush(stdout);
                exit(-1);
	}
	return(p);
}

static Kmer KmerAnd(Kmer kmer1, Kmer kmer2)
{
	kmer1.high1 &= kmer2.high1;
	kmer1.low1 &= kmer2.low1;
	kmer1.high2 &= kmer2.high2;
	kmer1.low2 &= kmer2.low2;
	return kmer1;
}

/*
static bool KmerEqualM(Kmer kmer1,Kmer kmer2)
{
	if(kmer1.high==kmer2.high&&kmer1.low==kmer2.low)
		return 1;
	else return 0;
}*/

static bool KmerEqualM(Kmer kmer1,Kmer kmer2)
{
	if(kmer1.high1==kmer2.high1&&kmer1.low1==kmer2.low1&&kmer1.high2==kmer2.high2&&kmer1.low2==kmer2.low2)
		return 1;
	else return 0;
}

static int searchKmerOnCtg(Kmer currW,Kmer *kmerDest,unsigned num)
{
	unsigned i;
	for(i=0;i<num;i++){
		if(KmerEqualM(currW,kmerDest[i])){
			return i;
		}
	}
	return -1;
}

static Kmer createFilter(int overlaplen)
{
	Kmer word;
	word.high1 = word.low1 = word.high2 = word.low2 = 0;
	if(2*overlaplen<64)
		word.low2 = (((ubyte8) 1)<<(2*overlaplen)) - 1;
	if(2*overlaplen>=64 && 2*overlaplen<128){
		word.low2 = ~word.low2;
		word.high2 = (((ubyte8) 1)<<(2*overlaplen-64)) - 1;
	}
	if(2*overlaplen>=128 && 2*overlaplen<192){
		word.low2 = ~word.low2;
		word.high2 = ~word.high2;
		word.low1 = (((ubyte8) 1)<<(2*overlaplen-128)) - 1;
	}
	if(2*overlaplen>=192 && 2*overlaplen<256){ 
		word.low2 = ~word.low2;
		word.high2 = ~word.high2;
		word.low1 = ~word.low1;
		word.high1 = (((ubyte8) 1)<<(2*overlaplen-192)) - 1; 
	}
	return word;
}

static void printKmerSeq(Kmer kmer,int overlaplen,int a)
{
	int i,bit1,bit2,bit3,bit4;
	char ch;
	char kmerSeq[128];
	if(overlaplen < 32) {bit4 = overlaplen; bit3 =0;bit2 =0;bit1 =0;}
	if(overlaplen >= 32 && overlaplen < 64) 
{bit4 = 32; bit3 =overlaplen-32;bit2 =0;bit1 =0;}
	if(overlaplen >= 64 && overlaplen < 96) 
{bit4 = 32; bit3 =32;bit2 =overlaplen-64;bit1 =0;}
	if(overlaplen >= 96 && overlaplen < 128) 
{bit4 = 32; bit3 =32;bit2 =32;bit1 =overlaplen-96;} 

	for(i=bit1-1;i>=0;i--){
		ch = kmer.high1&0x3;
		kmer.high1 >>= 2;
		kmerSeq[i] = ch;
	}
	for(i=bit2-1;i>=0;i--){
		ch = kmer.low1&0x3;
		kmer.low1 >>= 2;
		kmerSeq[i+bit1] = ch;
	}
	for(i=bit3-1;i>=0;i--){
		ch = kmer.high2&0x3;
		kmer.high2 >>= 2;
		kmerSeq[i+bit1+bit2] = ch;
	}
	for(i=bit4-1;i>=0;i--){
		ch = kmer.low2&0x3;
		kmer.low2 >>= 2;
		kmerSeq[i+bit1+bit2+bit3] = ch;
	}
	for(i=0;i<overlaplen;i++){
		//fprintf(fp,"%c",int2base((int)kmerSeq[i]));
		if (a)
		{
			cerr<<int2base((int)kmerSeq[i]);
		}else{
			cout<<int2base((int)kmerSeq[i]);
		}
	}
}

ubyte4 startPosition(){
	kset.iter_ptr=0;
	Kmer_t *rs;
	ubyte4 c=0;
	while (kset.iter_ptr < kset.size)
	{
		if (!is_kmer_entity_null(kset.flags, kset.iter_ptr))
		{
			rs=kset.array + kset.iter_ptr;
			ubyte4 num=get_kmer_num(*rs);
			set_kmer_start(*rs,c);
			c+=num;
		}
		kset.iter_ptr++;
	}
	return c;
}

void resetStart(){
	kset.iter_ptr=0;
	Kmer_t *rs;
	while (kset.iter_ptr < kset.size)
	{
		if (!is_kmer_entity_null(kset.flags, kset.iter_ptr))
		{
			rs=kset.array + kset.iter_ptr;
			ubyte4 num=get_kmer_num(*rs);
			ubyte4 start=get_kmer_start(*rs);
			set_kmer_start(*rs,start+1-num);
		}
		kset.iter_ptr++;
	}
}

void printR2K(vector<ubyte4>& readTkmer){
	Kmer_t *rs;
	kset.iter_ptr=0;
	while (kset.iter_ptr < kset.size)
	{
		if (!is_kmer_entity_null(kset.flags, kset.iter_ptr))
		{
			rs=kset.array + kset.iter_ptr;
			Kmer seq=get_kmer_seq(*rs);
			printKmerSeq(seq,args.overlap,1);
			ubyte4 start=get_kmer_start(*rs);
			ubyte4 num=get_kmer_num(*rs);
			cerr<<"\t"<<start<<"\t"<<num;
			for (int i=start; i!=start+num; ++i)
			{
				cerr<<"\t"<<readTkmer[i];
			}
			cerr<<endl;
		}
		kset.iter_ptr++;
	}
}

void printHash()
{
	Kmer_t *rs;
	kset.iter_ptr=0;
//	cerr<<"kset.size= "<<kset.size<<endl;
//	int c=0;
	while (kset.iter_ptr < kset.size){
		if (!is_kmer_entity_null(kset.flags, kset.iter_ptr))
		{
			rs=kset.array + kset.iter_ptr;
			//	__uint128_t temp;
			//	temp=kset.Kmer2int128( get_kmer_seq(*rs) );
//			c++;
			//cerr<<rs->seq.high<<rs->seq.low<<"\t";
			for (unsigned i=0; i<4; i++)
			{
				if (i<2)
				{
					cerr<<get_kmer_left_cov1(*rs,i)<<"\t";
				}else{
					cerr<<get_kmer_left_cov2(*rs,i-2)<<"\t";
				}
				
			}
			printKmerSeq(rs->seq,args.overlap,1);
			cerr<<"\t";
			for (unsigned i=0; i<4; i++)
			{
				if (i<2)
				{
					cerr<<get_kmer_right_cov1(*rs,i)<<"\t";
				}else{
					cerr<<get_kmer_right_cov2(*rs,i-2)<<"\t";
				}
			}
			cerr<<endl;
			// os1<<rs->checked<<"\t"<<rs->linear<<"\t"<<rs->refLinear<<"\t"<<rs->deleted<<"\t"<<rs->index<<"\t"<<rs->flag<<"\t"<<rs->single<<endl;
		}
		kset.iter_ptr++;
	}
//	cerr<<"Out kmer number is "<<c<<endl;
}

static bool KmerSmaller(Kmer kmer1,Kmer kmer2)
{
	if(kmer1.high1<kmer2.high1)
		return 1;
	else if(kmer1.high1==kmer2.high1){
		if(kmer1.low1<kmer2.low1)	
			return 1;
		else if(kmer1.low1==kmer2.low1){
			if(kmer1.high2<kmer2.high2)	
				return 1;
			else if(kmer1.high2==kmer2.high2){
				if(kmer1.low2<kmer2.low2)	
					return 1;
				else return 0; 
			}else return 0;
		}else return 0;
	}else return 0;
}

static Kmer KmerLeftBitMoveBy2(Kmer word)
{ 
	ubyte8 temp3 = word.low2>>62;
	ubyte8 temp2 = word.high2>>62;
	ubyte8 temp1 = word.low1>>62;
	word.high1 <<= 2;
	word.high1 |= temp1;
	word.low1 <<= 2;
	word.low1 |= temp2;
	word.high2 <<= 2;
	word.high2 |= temp3;
	word.low2 <<= 2;
	return word;
}

static Kmer KmerRightBitMoveBy2(Kmer word)
{
	ubyte8 temp1 = (word.high1&0x3)<<62;
	ubyte8 temp2 = (word.low1&0x3)<<62;
	ubyte8 temp3 = (word.high2&0x3)<<62;
	word.high1 >>= 2;
	word.low1 >>= 2;
	word.high2 >>= 2;  
	word.low2 >>= 2;  
	word.low1 |= temp1; 
	word.high2 |= temp2;
	word.low2 |= temp3; 
	return word;  
}

static Kmer nextKmer(Kmer prev,char ch)
{
	Kmer word = KmerLeftBitMoveBy2(prev);
	word = KmerAnd(word,WORDFILTER);
	word.low2 |= ch;
	return word;
}

static Kmer prevKmer(Kmer next,char ch)
{
	Kmer word = KmerRightBitMoveBy2(next);
	if(2*(args.overlap-1)<64)
		word.low2 |= (((ubyte8)ch) << 2*(args.overlap-1));
	if(2*(args.overlap-1)>=64 && 2*(args.overlap-1)<128)
		word.high2 |= ((ubyte8)ch) << (2*(args.overlap-1)-64);
	if(2*(args.overlap-1)>=128 && 2*(args.overlap-1)<192)
		word.low1 |= ((ubyte8)ch) << (2*(args.overlap-1)-128);
	if(2*(args.overlap-1)>=192 && 2*(args.overlap-1)<256)
		word.high1 |= ((ubyte8)ch) << (2*(args.overlap-1)-192);
	return word;
}
/*
static Kmer prevKmer(Kmer next,char ch)
{
	Kmer word = KmerRightBitMoveBy2(next);
	if(2*(args.overlap-1)<64)
		word.low |= (((ubyte8)ch) << 2*(args.overlap-1));
	else
		word.high |= ((ubyte8)ch) << (2*(args.overlap-1)-64);

	return word;
}

static bool KmerLarger(Kmer kmer1,Kmer kmer2)
{
	if(kmer1.high>kmer2.high)
		return 1;
	else if(kmer1.high==kmer2.high){
		if(kmer1.low>kmer2.low)	
			return 1;
		else 
			return 0;
	}else
		return 0;
}
*/
static bool KmerLarger(Kmer kmer1,Kmer kmer2)
{
	if(kmer1.high1>kmer2.high1)
		return 1;
	else if(kmer1.high1==kmer2.high1){
		if(kmer1.low1>kmer2.low1)	
			return 1;
		else if(kmer1.low1==kmer2.low1){
			if(kmer1.high2>kmer2.high2)	
				return 1;
			else if(kmer1.high2==kmer2.high2){
				if(kmer1.low2>kmer2.low2)	
					return 1;
				else return 0; 
			}else return 0;
		}else return 0;
	}else return 0;
}

static char lastCharInKmer(Kmer kmer)
{
	return (char)(kmer.low2&0x3);
}

static char firstCharInKmer (Kmer kmer)
{
	switch (args.overlap)
	{
	case 1 ... 32:
		kmer.low2 >>= 2 * (args.overlap - 1);
		return kmer.low2;	// & 3;
	case 33 ... 64:
		kmer.high2 >>= 2 * (args.overlap - 1) - 64;
		return kmer.high2;	// & 3;
	case 65 ... 96:
		kmer.low1 >>= 2 * (args.overlap - 1) - 128;
		return kmer.low1;
	case 97 ... 128:
		kmer.high1 >>= 2 * (args.overlap - 1) - 192;
		return kmer.high1;
	}
}

static Kmer KmerRightBitMove(Kmer word,int dis)
{
	if(dis<64){
		ubyte8 mask = (((ubyte8) 1) << dis) - 1;
		ubyte8 temp1 = (word.high1&mask)<<(64-dis);
		ubyte8 temp2 = (word.low1&mask)<<(64-dis);
		ubyte8 temp3 = (word.high2&mask)<<(64-dis);
		word.high1 >>= dis;
		word.low1 >>= dis;
		word.high2 >>= dis;
		word.low2 >>= dis;
		word.low1 |= temp1;
		word.high2 |= temp2;
		word.low2 |= temp3;
		return word;
	}
	if(dis>=64 && dis<128){
		ubyte8 mask = (((ubyte8) 1) << (dis-64)) - 1;
		ubyte8 temp1 = (word.high1&mask)<<(128-dis);
		ubyte8 temp2 = (word.low1&mask)<<(128-dis);
		word.high1 >>= (dis-64);
		word.low1 >>= (dis-64);
		word.high2 >>= (dis-64);
		word.high2 |= temp2;
		word.low2 = word.high2;
		word.low1 |= temp1;
		word.high2 = word.low1;
		word.low1 = word.high1;
		word.high1 = 0;
		return word;
	}
	if(dis>=128 && dis<192){
		ubyte8 mask = (((ubyte8) 1) << (dis-128)) - 1;
		ubyte8 temp1 = (word.high1&mask)<<(192-dis);
		word.high1 >>= (dis-128);
		word.low1 >>= (dis-128);
		word.low1 |= temp1;
		word.low2 = word.low1;
		word.high2 = word.high1;
		word.low1 =0;
		word.high1 = 0;
		return word;
	}
	if(dis>=192 && dis<256){
		word.high1 >>= (dis-192);
		word.low2 = word.high1;
		word.high2 = 0;
		word.low1 =0;
		word.high1 = 0;
		return word;
	}
}

static Kmer fastReverseComp(Kmer seq, char seq_size){
	seq.low2 ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.low2 = ((seq.low2 & 0x3333333333333333LLU)<< 2) | ((seq.low2 & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
	seq.low2 = ((seq.low2 & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seq.low2 & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
	seq.low2 = ((seq.low2 & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq.low2 & 0xFF00FF00FF00FF00LLU)>> 8);
	seq.low2 = ((seq.low2 & 0x0000FFFF0000FFFFLLU)<<16) | ((seq.low2 & 0xFFFF0000FFFF0000LLU)>>16);
	seq.low2 = ((seq.low2 & 0x00000000FFFFFFFFLLU)<<32) | ((seq.low2 & 0xFFFFFFFF00000000LLU)>>32);
	if(seq_size<32){
		seq.low2 >>= (64 - (seq_size<<1));
		return seq;
	}
	seq.high2 ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.high2 = ((seq.high2 & 0x3333333333333333LLU)<< 2) | ((seq.high2 & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
	seq.high2 = ((seq.high2 & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seq.high2 & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
	seq.high2 = ((seq.high2 & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq.high2 & 0xFF00FF00FF00FF00LLU)>> 8);
	seq.high2 = ((seq.high2 & 0x0000FFFF0000FFFFLLU)<<16) | ((seq.high2 & 0xFFFF0000FFFF0000LLU)>>16);
	seq.high2 = ((seq.high2 & 0x00000000FFFFFFFFLLU)<<32) | ((seq.high2 & 0xFFFFFFFF00000000LLU)>>32);
	ubyte8 temp;
	if(seq_size<64){
		temp = seq.high2;
		seq.high2 = seq.low2;
		seq.low2 = temp;
		seq = KmerRightBitMove(seq,128-(seq_size<<1));
		return seq;
	}
	seq.low1 ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.low1 = ((seq.low1 & 0x3333333333333333LLU)<< 2) | ((seq.low1 & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
	seq.low1 = ((seq.low1 & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seq.low1 & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
	seq.low1 = ((seq.low1 & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq.low1 & 0xFF00FF00FF00FF00LLU)>> 8);
	seq.low1 = ((seq.low1 & 0x0000FFFF0000FFFFLLU)<<16) | ((seq.low1 & 0xFFFF0000FFFF0000LLU)>>16);
	seq.low1 = ((seq.low1 & 0x00000000FFFFFFFFLLU)<<32) | ((seq.low1 & 0xFFFFFFFF00000000LLU)>>32);
	if(seq_size<96){
		temp = seq.low2;
		seq.low2 = seq.low1;
		seq.low1 = temp;
		seq = KmerRightBitMove(seq,192-(seq_size<<1));
		return seq;
	}
	seq.high1 ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.high1 = ((seq.high1 & 0x3333333333333333LLU)<< 2) | ((seq.high1 & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
	seq.high1 = ((seq.high1 & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seq.high1 & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
	seq.high1 = ((seq.high1 & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq.high1 & 0xFF00FF00FF00FF00LLU)>> 8);
	seq.high1 = ((seq.high1 & 0x0000FFFF0000FFFFLLU)<<16) | ((seq.high1 & 0xFFFF0000FFFF0000LLU)>>16);
	seq.high1 = ((seq.high1 & 0x00000000FFFFFFFFLLU)<<32) | ((seq.high1 & 0xFFFFFFFF00000000LLU)>>32);
	ubyte8 temp_t;
	temp = seq.high2;
	seq.high2 = seq.low1;
	seq.low1 = temp;
	temp_t = seq.high1;
	seq.high1 = seq.low2;
	seq.low2 = temp_t;
	seq = KmerRightBitMove(seq,256-(seq_size<<1));
	return seq;
}

static Kmer reverseComplement(Kmer word,int overlap)
{
	return fastReverseComp(word,overlap);
}

static void reverseComplementSeq(char *seq, int len,char *bal_seq)
{
	int i,index=0;
	
	if(len<1)
		return;

	for(i=len-1;i>=0;i--)
		bal_seq[index++] = int_comp(seq[i]);

	return;
}

static void chopKmer4End (Kmer *kmerCtg, int lenCtg, char *src_seq)
{
	int index, j;
	Kmer word;

	word.high1 = word.low1 = word.high2 = word.low2 = 0;

	for (index = 0; index < args.overlap; index++)
	{
		word = KmerLeftBitMoveBy2 (word);
		word.low2 |= src_seq[index];
	}

	index = 0;
	kmerCtg[index++] = word;

	for (j = 1; j <= lenCtg - args.overlap; j++)
	{
		word = nextKmer (word, src_seq[j - 1 + args.overlap]);
		kmerCtg[index++] = word;
	}
}

static void chopKmer4read(char *src_seq,unsigned len_seq,Kmer *kmerBuffer,uchar *prevcBuffer,uchar *nextcBuffer)
{
	int j,bal_j;
	Kmer word,bal_word;
	int index;
	char InvalidCh=4;
	char bal_seq[len_seq];

	if(len_seq<args.overlap+1){
		return;
	}
	word.high1 = word.low1 = word.high2 = word.low2 = 0;
	for (index = 0;index<args.overlap;index++){
		word = KmerLeftBitMoveBy2(word);
		word.low2 |= src_seq[index];
	}
	reverseComplementSeq(src_seq, len_seq,bal_seq);

		// complementary node
	bal_word = reverseComplement(word,args.overlap);
	bal_j = len_seq-0-args.overlap;  //  0;     
	index = 0;
	if(KmerSmaller(word,bal_word)){
		kmerBuffer[index] = word;
		prevcBuffer[index] = InvalidCh;
		nextcBuffer[index++] = src_seq[0+args.overlap];
	}else{
		kmerBuffer[index] = bal_word;
		prevcBuffer[index] = bal_seq[bal_j-1];
		nextcBuffer[index++] = InvalidCh;
	}
//	cerr<<"0\t"<<word.high<<word.low<<"\t"<<bal_word.high<<bal_word.low<<endl;
	for(j = 1; j <= len_seq - args.overlap; j ++)	{
		word = nextKmer(word,src_seq[j-1+args.overlap]); 
		bal_j = len_seq-j-args.overlap; //  j; 
		bal_word = prevKmer(bal_word,bal_seq[bal_j]);
		
		if(KmerSmaller(word,bal_word)){
			kmerBuffer[index] = word;
			prevcBuffer[index] = src_seq[j-1];
			if(j<len_seq - args.overlap)
				nextcBuffer[index++] = src_seq[j+args.overlap];
			else
				nextcBuffer[index++] = InvalidCh;
				//printf("%dth: %p with %p\n",kmer_c-1,word,hashBanBuffer[kmer_c-1]);
		}else{
				// complementary node
			kmerBuffer[index] = bal_word;
			if(bal_j>0)
				prevcBuffer[index] = bal_seq[bal_j-1];
			else
				prevcBuffer[index] = InvalidCh;
			nextcBuffer[index++] = bal_seq[bal_j+args.overlap];
				//printf("%dth: %p with %p\n",kmer_c-1,bal_word,hashBanBuffer[kmer_c-1]);
		}
//		cerr<<j<<"\t"<<word.high<<word.low<<"\t"<<bal_word.high<<bal_word.low<<endl;
	}
}

static void codeDBGraph(string& s,char* src){
	string::iterator its=s.begin();
	unsigned c=0;
	for (; its!=s.end(); ++its)
	{
		if (*its>='A' && *its<='Z')
		{
			src[c++]=base2int(*its);
		}else if (*its>='a' && *its<='z')
		{
			src[c++]=base2int(*its-'a'+'A');
		}else if (*its=='.')
		{
			src[c++]=base2int('A');
		}
	}
}

static void putKmer2DBGraph(Kmer *kmerBuffer,uchar *prevcBuffer,uchar *nextcBuffer,unsigned len){
	Kmer_t *node;
    for (unsigned j=0; j != len-args.overlap+1; ++j)
    {
//		cout<<kmerBuffer[j].high1<<kmerBuffer[j].low1<<kmerBuffer[j].high2<<kmerBuffer[j].low2<<"\t";
//		printKmerSeq(kmerBuffer[j],args.overlap,0);
//		cout<<"\t"<<int(prevcBuffer[j])<<"\t"<<int(nextcBuffer[j])<<endl;;
		kset.put_kmerset(kmerBuffer[j], prevcBuffer[j], nextcBuffer[j], &node);
	}
}

static void putReadIndex2Vector(Kmer *kmerBuffer,unsigned len,ubyte4 index){
	Kmer_t *node;
	for (unsigned j=0; j!=len-args.overlap+1; ++j)
	{
		bool found = kset.search_kmerset(kmerBuffer[j],&node);
		if (!found)
		{
			cerr<<"Err: can't find kmer " <<endl;
		}else{
			ubyte4 s=node->start; //node->start should points at the last plug-in readID index. ct
            if(readTkmer[s]==0){
                readTkmer[s] = index;//if nothing has been plug into, put the index into readTkmer
            }
            //check the readID index is repetitive. ct
            //if repetitive, not put into readTkmer and kmer.num subtract by 1. ct
            else if(readTkmer[s]==index){
                ubyte4 num = node->num;
                set_kmer_num(*node, num-1); 
            }else{    //otherwise, add into readTkmer, and move node->start down by 1
                readTkmer[s+1]=index;
                set_kmer_start(*node,s+1);
            }
			/*ubyte4 s=node->start;
			readTkmer[s]=index;
			set_kmer_start(*node,s+1);
			*/
		}
	}
}

static void chopRead2Kmer(vector<string>& temp,char flag,ubyte4 index){
	Kmer kmerBuffer[1024];
	uchar prevBuffer[1024];
	uchar nextBuffer[1024];
//	cout<<"chopRead2Kmer temp.size()="<<temp.size()<<endl;
	for (int it=0; it!=temp.size(); ++it)
	{
		
		unsigned len=temp[it].size();
		char src[len];
		codeDBGraph(temp[it],src);
		/*for (; its!=temp[it].end(); ++its)
		{
			if (*its>='A' && *its<='Z')
			{
				src[c++]=base2int(*its);
			}else if (*its>='a' && *its<='z')
			{
				src[c++]=base2int(*its-'a'+'A');
			}else if (*its=='.')
			{
				src[c++]=base2int('A');
			}
		}
		for (unsigned i=0; i!=len; ++i)
		{
			cout<<int(src[i]);
		}
		cout<<endl;*/
		chopKmer4read(src,len,kmerBuffer,prevBuffer,nextBuffer);
//		cout<<it<<"\tafter chopKmer4read"<<endl;
		if (flag)
		{
			putKmer2DBGraph(kmerBuffer,prevBuffer,nextBuffer,len);
		}else{
			putReadIndex2Vector(kmerBuffer,len,index+it+1);
		}
//		cout<<it<<"\tafter putKmer2DBGraph"<<endl;
	}
}//void chopKmer4read(char *src_seq,unsigned len_seq,Kmer *kmerBuffer,uchar *prevcBuffer,uchar *nextcBuffer)

static void readFasta(ifstream& handle,char flag){
	int a=0;int b=0;
	string line;
	vector<string> temp;
	while (getline(handle,line))
	{
		getline(handle,line);
		++a;
		temp.push_back(line);
		if (a==1000000)
		{
			/*for (unsigned i=0; i!=temp.size(); ++i)
			{
				cout<<temp[i]<<endl;
			}*/
			chopRead2Kmer(temp,flag,b);
			temp.clear();
			b+=a;
			a=0;
		}
	}
	if (!temp.empty())
	{
		/*for (unsigned i=0; i!=temp.size(); ++i)
		{
			cout<<temp[i]<<endl;
		}*/
		chopRead2Kmer(temp,flag,b);
	}
//	handle.close();
}

static void readFastq(ifstream& handle,char flag){
	string line;
	vector<string> temp;
	int a=0;
	int b=0;
	while (getline(handle,line))
	{
		getline(handle,line);
		++a;
//		cout<<a<<"\t"<<line<<endl;
		temp.push_back(line);
		if (a==1000000) //1000000
		{
			/*for (unsigned i=0; i!=temp.size(); ++i)
			{
				cout<<temp[i]<<endl;
			}*/
//			++b;
			chopRead2Kmer(temp,flag,b);
			temp.clear();
			b+=a;
			a=0;
//			cout<<"b = "<<b<<endl;
			//break;
			//exit(1);
		}
		getline(handle,line);
		if (line[0]!='+')
		{
			cout<<"Error data formate!\n";
			exit(1);
		}
		getline(handle,line);
	}
//	cout<<"readFastq 3\n";
	if (!temp.empty())
	{
		chopRead2Kmer(temp,flag,b);
	}
//	handle.close();
}

static void inputEndsFile(ifstream& handle,vector<string>& end1,vector<string>& end2){
	/*ends file formate must be as below,
	>id
	end1
	end2*/
	string line;
	while (getline(handle,line))
	{
		if (line[0]!='>')
		{
			cout<<"wrong ends formate!"<<endl;
			exit(-1);
		}
//		cout<<line.size()<<"\t"<<line<<endl;
		getline(handle,line);
		end1.push_back(line);
//		cout<<line.size()<<"\t"<<line<<endl;
		getline(handle,line);
		end2.push_back(line);
//		cout<<line<<endl;
	}
	if (end1.size()!=end2.size())
	{
		cout<<"error while input ends file"<<endl;
		exit(-1);
	}
}
//##########
static unsigned maxRightCov(Kmer_t *node ){
	char ch;
	unsigned links,links1=0;
	for(ch=0;ch<4;ch++){
		if (ch<2)
		{
			links = get_kmer_right_cov1(*node,ch);
		}else{
			links = get_kmer_right_cov2(*node,ch-2);
		}
		
		if (links1<links)
		{
			links1=links;
		}
	}
	return links1;
}

static unsigned maxLeftCov(Kmer_t *node ){
	unsigned ch;
	unsigned links,links1=0;
	for(ch=0;ch<4;ch++){
		if (ch<2)
		{
			links = get_kmer_left_cov1(*node,ch);
		}else{
			links = get_kmer_left_cov2(*node,ch-2);
		}
		
		if (links1<links)
		{
			links1=links;
		}
	}
	return links1;
}
/*
static int clipTipFromNode(Kmer_t *node1,int cut_len)
{
	uchar ret=0,in_num,out_num,link;
	int sum,count,total=0;
	Kmer_t *out_node,*pre_node;
	Kmer tempKmer,pre_word,word,bal_word;
	char ch1,ch;
	bool smaller,found;
//	int setPicker;
	unsigned int max_links,singleCvg;

	if(node1->linear || node1->deleted)
		return ret;
	if(node1->checked)
		return ret;
	in_num = kset.count_branch2prev(node1);
	out_num = kset.count_branch2next(node1);
	if(in_num==0&&out_num==1){
		pre_word = node1->seq;
		for(ch1=0;ch1<4;ch1++){
			link = get_kmer_right_cov(*node1,ch1);
			if(link)
				break;
		}
		word = nextKmer(pre_word,ch1);
	}else if(in_num==1&&out_num==0){
		pre_word = fastReverseComp(node1->seq,args.overlap);
		for(ch1=0;ch1<4;ch1++){
			link = get_kmer_left_cov(*node1,ch1);
			if(link)
				break;
		}
		word = nextKmer(pre_word,int_comp(ch1));  
	}else
		return ret;

	count = 1;
	bal_word = fastReverseComp(word,args.overlap);
	if(KmerLarger(word,bal_word) ){ //word>bal_word
		tempKmer = bal_word;
		bal_word = word;
		word = tempKmer;
		smaller = 0;
	}else
		smaller = 1;
	
	pre_node=node1;
	found = kset.search_kmerset( word, &out_node);
	if(!found){
		printf("kmer %llx not found, node1 %llx\n",word,node1->seq);
		exit(1);
	}
	while(out_node->linear){
		count++;
		if(out_node->checked)
			return ret;
		if(count>cut_len)
			return ret;
		if(smaller){
			pre_word = word;
			for(ch=0;ch<4;ch++){
				link = get_kmer_right_cov(*out_node,ch);
				if(link)
					break;
			}

			word = nextKmer(pre_word,ch);
			bal_word = fastReverseComp(word,args.overlap);
			if(KmerLarger(word,bal_word)){ //word>bal_word
				tempKmer = bal_word;
				bal_word = word;
				word = tempKmer;
				smaller = 0;
			}else
				smaller = 1;
			
			pre_node=out_node;
			found = kset.search_kmerset( word, &out_node);
			if(!found){
				printf("kmer %llx not found, node1 %llx\n",word,node1->seq);
				printf("pre_word %llx with %d(smaller)\n",pre_word,ch);
				exit(1);
			}
		}
		else{
			pre_word = bal_word;
			for(ch=0;ch<4;ch++){
				link = get_kmer_left_cov(*out_node,ch);
				if(link)
					break;
			}

			word = nextKmer(pre_word,int_comp(ch));
			bal_word = fastReverseComp(word,args.overlap);
			if(KmerLarger(word,bal_word)){
				tempKmer = bal_word;
				bal_word = word;
				word = tempKmer;
				smaller = 0;
			}else
				smaller = 1;
			
			pre_node=out_node;
			found = kset.search_kmerset( word, &out_node);
			if(!found){
				printf("kmer %llx not found, node1 %llx\n",word,node1->seq);
				printf("pre_word %llx with %d(larger)\n",fastReverseComp(pre_word,args.overlap),int_comp(ch));
				exit(1);
			}
		}
	}

	if((sum = kset.count_branch2next(out_node) + kset.count_branch2prev(out_node)) ==1){
		tip_c++;
		node1->deleted = 1;
		out_node->deleted = 1;
		return 1;
	}else
	{
		ch = firstCharInKmer(pre_word);
		// make sure this tip doesn't provide most links to out_node
		max_links = 0;
		for(ch1=0;ch1<4;ch1++){
			if(smaller){
				singleCvg = get_kmer_left_cov(*out_node,ch1);	
				if(singleCvg>max_links)
					max_links = singleCvg;
			}else{
				singleCvg = get_kmer_right_cov(*out_node,ch1);	
				if(singleCvg>max_links)
					max_links = singleCvg;
			}
		}
		
		if(smaller && (get_kmer_left_cov(*out_node,ch)<max_links)){
		//if(smaller ){
			tip_c++;
			node1->deleted = 1;
			pre_node->deleted = 1;
			kset.dislink2prevUncertain(out_node,ch,smaller);
			if(kset.count_branch2prev(out_node)==1&&
					kset.count_branch2next(out_node)==1)
				out_node->linear = 1;
			return 1;
		}
		if(!smaller && (get_kmer_right_cov(*out_node,int_comp(ch))<max_links)){
		//if(!smaller){
			tip_c++;
			node1->deleted = 1;
			pre_node->deleted = 1;
			kset.dislink2prevUncertain(out_node,ch,smaller);
			if(kset.count_branch2prev(out_node)==1&&
					kset.count_branch2next(out_node)==1)
				out_node->linear = 1;
			return 1;
		}
	}
	return 0;
}

void removeTips()
{
	int i,flag=1,cut_len_tip;
	Kmer_t *rs;
		
	//count_ends(hash_table);
	//cut_len_tip = 2*overlaplen >= maxReadLen4all-overlaplen+1 ? 2*overlaplen : maxReadLen4all-overlaplen+1; 
	cut_len_tip = 2*args.overlap;
	//printf("Start to remove tips of single frequency kmers short than %d\n",cut_len_tip);
	tip_c = 0;
	while(flag){			//
		flag = 0;
		kset.iter_ptr = 0;
		while(kset.iter_ptr < kset.size){
			if(!is_kmer_entity_null(kset.flags, kset.iter_ptr)){
				rs = kset.array + kset.iter_ptr;
				flag += clipTipFromNode(rs,cut_len_tip);
			}
			kset.iter_ptr ++;
		}
	}
	cout<<"Remove "<<tip_c<<" tips"<<endl;
}
*/
int abundanceSort(Kmer_t *node,char *array,int cutlink,char flag){
	char ch;
	unsigned a[4]={0,0,0,0};
	int c=0;
	unsigned links;
	if (flag)
	{
		for (ch=0; ch<4; ch++)
		{
			if (ch<2)
			{
				links = get_kmer_right_cov1(*node,ch);
			}else{
				links = get_kmer_right_cov2(*node,ch-2);
			}
			a[ch]=links;
			if (links>=cutlink)
			{
				++c;
			}
		}
	}else{
		for (ch=0; ch<4; ch++)
		{
			if (ch<2)
			{
				links = get_kmer_left_cov1(*node,ch);
			}else{
				links = get_kmer_left_cov2(*node,ch-2);
			}
			a[ch]=links;
			if (links>=cutlink)
			{
				++c;
			}
		}
	}
/*	cerr<<"c\t"<<c<<"\t";
	for (ch=0; ch!=4; ch++)
	{
		cerr<<a[ch]<<"\t";
	}
	cerr<<endl;
*/
	char b1,b2,b3,b4,b5,b6;
	b1=a[0]>a[1]?0:1;
	b2=a[2]>a[3]?2:3;
	b3=1-b1;
	b4=5-b2;
	if (a[b1]>a[b2])
	{
		array[0]=b1;
		b5=b2;
	}else{
		array[0]=b2;
		b5=b1;
	}
	if (a[b3]<a[b4])
	{
		array[3]=b3;
		b6=b4;
	}else{
		array[3]=b4;
		b6=b3;
	}
	if (a[b5]>a[b6])
	{
		array[1]=b5;
		array[2]=b6;
	}else{
		array[1]=b6;
		array[2]=b5;
	}
//	cerr<<"b\t"<<(int)b1<<"\t"<<(int)b2<<"\t"<<(int)b3<<"\t"<<(int)b4<<"\t"<<(int)b5<<"\t"<<(int)b6<<endl;
	return c;
}

char resoveBranch(Kmer_t **node,Kmer_t **soFarNode,vector<int>& branchArr,int steps,int cutSup){
//	cerr<<"resoveBranch 1 "<<endl;
	int i,j;
	Kmer_t *node1;
	ubyte4 start,num;
	i=steps-1; //add "-1"
	unordermap m1,m2,m3;
//	cerr<<"resoveBranch 2 step= "<<steps<<endl;
	for (; i>0; i--) //i mey be used to control whether is the first branch
	{
		if (branchArr[i]>1)
		{
			break;
		}
	}
	unordermap::iterator it;
	int freq=0;
	int k1=i;int k2=0;
	if (i)
	{
		if (steps-i>30)
		{
			i=steps-30;
		}
		k2=i;
		node1=soFarNode[i-1];//i-1 error
//		cerr<<"resoveBranch 3.1 "<<node1<<"\t"<<node<<endl;	

		start=get_kmer_start(*node1);
		num=get_kmer_num(*node1);
//		cerr<<"resoveBranch 3.2 start = "<<start<<"\tnum = "<<num<<endl;
		for (j=start; j!=start+num; j++)
		{
//			cerr<<"resoveBranch 3.3 "<<j<<"\t"<<readTkmer[j]<<endl;
			m1[readTkmer[j] ]++;
		}
		
//		cerr<<"resoveBranch 4 "<<endl;
		for (; i<steps; i++)//=steps?
		{
			node1=soFarNode[i];
			start=get_kmer_start(*node1);
			num=get_kmer_num(*node1);
			for (j=start; j!=start+num; j++)
			{
				m2[readTkmer[j] ]++;
			}
		}
//		cerr<<"resoveBranch 5 "<<endl;
		node1=soFarNode[steps];
		start=get_kmer_start(*node1);
		num=get_kmer_num(*node1);
		for (j=start; j!=start+num; j++)
		{
			m3[readTkmer[j] ]++;
		}
		for (it=m1.begin(); it!=m1.end(); ++it)
		{
			if (m2.find(it->first)!=m2.end())
			{
				if (m3.find(it->first)!=m3.end())
				{
					freq++;
				}	
			}
		}
		cerr<<"^\t"<<m1.size()<<"\t"<<m2.size()<<"\t"<<m3.size()<<"\t"<<freq<<"\t"<<k1<<"\t"<<k2<<"\t"<<steps<<"\t"<<cutSup<<endl;
	}else{
		if (steps-i>30)
		{
			//steps and steps-50
			node1=soFarNode[steps-30];
			start=get_kmer_start(*node1);
			num=get_kmer_num(*node1);
			for (j=start; j!=start+num; j++)
			{
				m1[readTkmer[j] ]++;
			}
			node1=soFarNode[steps];
			start=get_kmer_start(*node1);
			num=get_kmer_num(*node1);
			for (j=start; j!=start+num; j++)
			{
				m2[readTkmer[j] ]++;
			}
		}else{
			//steps and i
			node1=soFarNode[i];
			start=get_kmer_start(*node1);
			num=get_kmer_num(*node1);
			for (j=start; j!=start+num; j++)
			{
				m1[readTkmer[j] ]++;
			}
			node1=soFarNode[steps];
			start=get_kmer_start(*node1);
			num=get_kmer_num(*node1);
			for (j=start; j!=start+num; j++)
			{
				m2[readTkmer[j] ]++;
			}
		}
		for (it=m1.begin(); it!=m1.end(); ++it)
		{
			if (m2.find(it->first)!=m2.end())
			{
				freq++;	
			}
		}
		cerr<<"#\t"<<m1.size()<<"\t"<<m2.size()<<"\t"<<m3.size()<<"\t"<<freq<<"\t"<<i<<"\t"<<steps<<"\t"<<cutSup<<endl;
	}

//	cerr<<"resoveBranch 3 i="<<i<<endl;
//	i=i-1;
//	cerr<<"resoveBranch 6 "<<endl;
	
//	cerr<<"resoveBranch 7 "<<endl;
	//if (freq>args.readSup) //this value should be tested args.readSup
	if (freq>cutSup)
	{
		return 1;
	}else{
		return 0;
	}
//	cerr<<"resoveBranch 8 "<<endl;
}

static void e2eass(Kmer currW,Kmer *kmerDest,int steps,int *num_route,
			int *traceCounter,long long *soFarLinks,Kmer_t **soFarNode,vector<int>& routeLens,
			vector<int>& routeEndOnE2,vector<char>& soFarSeq,vector<vector<char> >& foundRoutes,
			vector<double> avgLinks,unsigned num,vector<int>& branchArr,int branches,vector<int>& traceSteps)
{
	(*traceCounter)++;
	if(*traceCounter>UPlimit)
		return;
	if(steps>args.max||*num_route>=MaxRouteNum)
		return;
	Kmer word = reverseComplement(currW,args.overlap);
	bool isSmaller = KmerSmaller(currW,word);
	char ch,ch1;
	unsigned links,links1;
	int index,i;

	if(isSmaller)
		word = currW;
	
	Kmer_t *node;
	printKmerSeq(currW,args.overlap,1);
	cerr<<"\t";
	printKmerSeq(word,args.overlap,1);
	cerr<<"\t"<<steps<<"\t"<<isSmaller;
	bool found = kset.search_kmerset(word,&node);
	if(!found){
		cout<<"Trace: can't find kmer\t";
		printKmerSeq(word,args.overlap,0);
		cout<<"\t";
		printKmerSeq(currW,args.overlap,0);
		cout<<"\t"<<steps<<endl;
		return;
	}
	soFarNode[steps] = node;
	branchArr[steps]=branches;
	int linkCounter,cutlink,cutSup,cutTemp;
	if(soFarLinks){
		linkCounter = *soFarLinks;
		cutTemp=(int)(linkCounter/(steps+1));
		cutlink=(int)(cutTemp*args.lowFreCut);
		cutSup=(int)(cutTemp*args.readSup); 
	}
	cerr<<"\t"<<linkCounter<<"\t"<<cutTemp<<endl;
//	cerr<<"e2eass node add= "<<node<<"\t"<<steps<<endl;
	if (branches>1)
	{
		char ff=resoveBranch(&node,soFarNode,branchArr,steps,cutSup);
		if (!ff)
		{
			return;
		}
	}
//	branchArr[steps]=branches;
//	if(soFarNode)
//		soFarNode[steps] = node;
	if(steps>0)
		soFarSeq[steps-1] = lastCharInKmer(currW);
	
	int end;
	if(steps>=args.min&&((end=searchKmerOnCtg(word,kmerDest,num))>=0)){
		cerr<<"####\t";
		printKmerSeq(currW,args.overlap,1);
		cerr<<"\t";
		printKmerSeq(word,args.overlap,1);
		cerr<<"\t";
		printKmerSeq(kmerDest[end],args.overlap,1);
		cerr<<endl;
		index = *num_route; 
		if(steps>0)
			avgLinks[index] = (double)linkCounter/steps;
		else
			avgLinks[index] = 0;
		
		routeLens[index] = steps;
//		cerr<<"%%%%% the index "<<index<<" end is "<<end<<"\t"<<num<<endl;
		routeEndOnE2[index]=end;
		traceSteps[index]=*traceCounter;//for test
		for(i=0;i<steps;i++)
			foundRoutes[index][i]=soFarSeq[i];
		if(i<args.max)
			foundRoutes[index][i] = 4;  //indicate the end of the sequence
		*num_route = ++index;
		return;
	}

	steps++;
	char array[4]={0,0,0,0};
	if(isSmaller){
		branches=abundanceSort(node,array,cutlink,1);

/*		cerr<<"Array: ";
		for(ch1=0;ch1<4;ch1++){
			cerr<<int(array[ch1])<<"\t";
		}

		for (ch1=0; ch1<4; ch1++)
		{
			ch=array[ch1];
			if (ch<2)
			{
				links = get_kmer_right_cov1(*node,ch);
			}else{
				links = get_kmer_right_cov2(*node,ch-2);
			}
			cerr<<(int)links<<"\t";
		}
		cerr<<endl;
*/
		for(ch1=0;ch1<4;ch1++){
			ch=array[ch1];
			if (ch<2)
			{
				links = get_kmer_right_cov1(*node,ch);
			}else{
				links = get_kmer_right_cov2(*node,ch-2);
			}
			if (links<cutlink || links<1)
			{
				continue;
			}
		
//			if(links1<args.lowFreCut)
//				continue;
			cerr<<"+\t"<<steps<<"\t"<<int(links)<<"\t"<<cutlink<<"\t"<<int(ch)<<"\t"<<int2base(ch)<<"\t"<<branches<<endl;
			if(soFarLinks)
				*soFarLinks = linkCounter + links;
			word = nextKmer(currW,ch);
			e2eass(word,kmerDest,steps,num_route,traceCounter,soFarLinks,soFarNode,routeLens,routeEndOnE2,
				soFarSeq,foundRoutes,avgLinks,num,branchArr,branches,traceSteps);
		}
	}else{
/*		cerr<<"links: ";
		for (ch1=0;ch1<4;ch1++ )
		{
			if (ch1<2)
			{
				links = get_kmer_left_cov1(*node,ch1);
			}else{
				links = get_kmer_left_cov2(*node,ch1-2);
			}
			cerr<<int(links)<<"\t";
		}
*/
		branches=abundanceSort(node,array,cutlink,0);
/*		cerr<<"Array: ";
		for(ch1=0;ch1<4;ch1++){
			cerr<<int(array[ch1])<<"\t";
		}
		for (ch1=0; ch1<4; ch1++)
		{
			ch=array[ch1];
			if (ch<2)
			{
				links = get_kmer_left_cov1(*node,ch);
			}else{
				links = get_kmer_left_cov2(*node,ch-2);
			}
			cerr<<(int)links<<"\t";
		}
		cerr<<endl;
*/		
		for(ch1=0;ch1<4;ch1++){
			ch=array[ch1];
			if (ch<2)
			{
				links = get_kmer_left_cov1(*node,ch);
			}else{
				links = get_kmer_left_cov2(*node,ch-2);
			}
			if (links<cutlink || links<1)
			{
				continue;
			}
		 
//			if(links1<args.lowFreCut)
//				continue;
			cerr<<"-\t"<<steps<<"\t"<<int(links)<<"\t"<<cutlink<<"\t"<<int(ch)<<"\t"<<int2base(ch)<<"\t"<<branches<<endl;
			if(soFarLinks)
				*soFarLinks = linkCounter + links;
			word = nextKmer(currW,int_comp(ch));
			e2eass(word,kmerDest,steps,num_route,traceCounter,soFarLinks,soFarNode,routeLens,routeEndOnE2,
				soFarSeq,foundRoutes,avgLinks,num,branchArr,branches,traceSteps);
		}
	}
}
//#############
/*
void e2eass(kmerBuffer1[i],kmerBuffer2){
}
static string convertCode(string& s){
	string str;
	for (indexT i=0; i!=s.size(); ++i)
	{
		str+=alph.encode[s[i] ];
	}
	return str;
}*/

inline long long covStartKmer(Kmer prev,Kmer curr){
	char c=lastCharInKmer(curr);
	Kmer word = reverseComplement(prev,args.overlap);
	bool isSmaller = KmerSmaller(prev,word);
	if (isSmaller)
	{
		word=prev;
	}
	Kmer_t *node;
//	printKmerSeq(currW,args.overlap,1);
//	cerr<<"\t";
//	printKmerSeq(word,args.overlap,1);
//	cerr<<"\t"<<steps<<"\t"<<isSmaller<<endl;
	bool found = kset.search_kmerset(word,&node);
	if (!found)
	{
		cerr<<"ERR:Can't find end sequence kmer!\n";
		exit(1);
	}
//	cerr<<"covStartKmer "<<int(c)<<endl; 
	unsigned links;
	if(isSmaller){
		/*
		cerr<<"111"<<endl;
		for(int ch1=0;ch1<4;ch1++){
			if (ch1<2)
			{
				links = get_kmer_right_cov1(*node,ch1);
			}else{
				links = get_kmer_right_cov2(*node,ch1-2);
			}
			cerr<<int(links)<<"\t";
		}
		cerr<<endl;
		*/
		if (c<2)
		{
				return get_kmer_right_cov1(*node,c);
		}else{
				return get_kmer_right_cov2(*node,c-2);
		}
	}else{
		/*
		cerr<<"222"<<endl;
		for(int ch1=0;ch1<4;ch1++){
			if (ch1<2)
			{
				links = get_kmer_left_cov1(*node,ch1);
			}else{
				links = get_kmer_left_cov2(*node,ch1-2);
			}
			cerr<<int(links)<<"\t";
		}
		cerr<<endl;
		*/
		c=int_comp(c);
		if (c<2)
		{
			return get_kmer_left_cov1(*node,c);
		}else{
			return get_kmer_left_cov2(*node,c-2);
		}
	}
}

static void ass1by1(vector<string>& end1,vector<string>& end2,int index){
//	cerr<<"ass1by1 1 "<<endl;
	unsigned e1,e2;
	e1=end1[index].size();e2=end2[index].size();
	char src1[e1];
	char src2[e2];
	
	codeDBGraph(end1[index],src1);
	codeDBGraph(end2[index],src2);
	
	Kmer kmerBuffer1[1024];
	Kmer kmerBuffer2[1024];
	uchar prevBuffer[1024];
	uchar nextBuffer[1024];
//	cout<<"ass1by1 2 end1 "<<endl;
	chopKmer4read(src1,e1,kmerBuffer1,prevBuffer,nextBuffer);
	/*cout<<"end1 "<<e1<<"\t"<<end1[index]<<endl;
	for (unsigned j=0; j != e1-args.overlap+1; ++j)
    {
		printKmerSeq(kmerBuffer1[j],args.overlap,0);
		cout<<"\t"<<int(prevBuffer[j])<<"\t"<<int(nextBuffer[j])<<endl;;
	}*/
	putKmer2DBGraph(kmerBuffer1,prevBuffer,nextBuffer,e1);
	chopKmer4End (kmerBuffer1, e1, src1);
//	cout<<"ass1by1 2.1 end2 "<<endl;
	chopKmer4read(src2,e2,kmerBuffer2,prevBuffer,nextBuffer);
	/*cout<<"end2 "<<e2<<"\t"<<end2[index]<<endl;
	for (unsigned j=0; j != e2-args.overlap+1; ++j)
    {
		printKmerSeq(kmerBuffer2[j],args.overlap,0);
		cout<<"\t"<<int(prevBuffer[j])<<"\t"<<int(nextBuffer[j])<<endl;;
	}*/
	putKmer2DBGraph(kmerBuffer2,prevBuffer,nextBuffer,e2);
//	cout<<"ass1by1 3 "<<endl;
	
//	vector<Kmer_t *> soFarNode(args.max+1);
	Kmer_t **soFarNode;
	soFarNode=(Kmer_t **)ckalloc((args.max+1)*sizeof(Kmer_t *));
	vector<int> routeLens(MaxRouteNum);
	vector<char> soFarSeq(args.max);
	vector<vector<char> > foundRoutes(MaxRouteNum,vector<char>(args.max));
	vector<double> avgLinks(MaxRouteNum);
	vector<int> routeEndOnE2(MaxRouteNum);
	vector<int> branchArr(args.max+1);
//	cerr<<"ass1by1 4 "<<endl;
	
//	cerr<<end1[index].size()<<"\t"<<args.overlap<<"\t"<<end1[index].size()-args.overlap<<endl;
	for (int i=e1-args.overlap; i>=1; i-=10)
	{
//		cerr<<"ass1by1 4.1 "<<endl;
		int num_route=0;
		int traceCounter=0;
		long long soFarLinks=0;
		vector<int> traceSteps(MaxRouteNum);//for test
		int steps=0;
//		cerr<<"ass1by1 5 "<<i<<"\t";
		cerr<<">the "<<index<<" ed sequence, start at "<<i<<endl;
		printKmerSeq(kmerBuffer1[i],args.overlap,1);
		cerr<<endl;
		soFarLinks=covStartKmer(kmerBuffer1[i-1],kmerBuffer1[i]);
		cerr<<"soFarLinks "<<soFarLinks<<endl;
		e2eass(kmerBuffer1[i],kmerBuffer2,steps,&num_route,&traceCounter,&soFarLinks, soFarNode,
				 routeLens,routeEndOnE2,soFarSeq,foundRoutes,avgLinks,e2-args.overlap+1,branchArr,0,traceSteps);
//		cerr<<"ass1by1 6 "<<routeLens.size()<<"\t"<<num_route<<endl;
		if (num_route>0)
		{
			//for (int j=0; j!=routeLens.size(); ++j)
			flagResult[index]=1;
			for (int j=0; j!=num_route; ++j)
			{
				string a,b,saCode;
				char c;
				for (int k=0; k!=args.max; ++k)
				{
					if (foundRoutes[j][k]>3)
					{
						break;
					}
					c=int2base(foundRoutes[j][k]);
					a+=c;
//					saCode+=alph.encode[c];
					//b+=int(foundRoutes[j][k]);
				}
//				cout<<">"<<index<<"\t"<<i<<"\t"<<routeEndOnE2[j]<<endl;
				string prevs,nexts;
				prevs=end1[index].substr(0,i+args.overlap);
				nexts=end2[index].substr(routeEndOnE2[j]+args.overlap,e2-routeEndOnE2[j]-args.overlap);
//				cout<<a<<endl;
//				cout<<prevs<<a<<nexts<<endl;
//				cout<<end1[index]<<endl;
//				cout<<end2[index]<<endl;
//				saCode=convertCode(prevs)+saCode;
//				saCode+=convertCode(nexts);
//				contig.push_back(saCode); //check whether need to add prevs and nexts or not
				CONTIG[index].push_back(prevs+a+nexts);
				traceCount[index].push_back(traceSteps[j]);
				indexT aLen=a.size();
//				cout<<aLen<<"\n"<<prevs<<a<<nexts<<endl;
				/*for (unsigned kk=0; kk!=saCode.size(); ++kk)
				{
					cout<<alph.decode[saCode[kk] ];
				}
				cout<<endl;*/
				//cout<<nexts<<endl;
				//cout<<b<<endl;
//				cout<<"ass1by1 7 "<<routeLens.size()<<"\t"<<j<<endl;
			}
//			cerr<<"ass1by1 7.1 "<<i<<endl;
			break;
		}
//		cerr<<"ass1by1 8 "<<i<<endl;
	}
	//exit(1);
//	cerr<<"ass1by1 9 "<<endl;
	
}

static void ass(int i,std::ostream& outs){
//	cout<<"begin at ass() "<<end1.size()<<"\t"<<end2.size()<<endl;
	
		cout<<"ass() begin at ends "<<i<<endl;
		vector<string> contig;
		ass1by1(end1,end2,i);
		stringstream stream;
		stream<<i;
		string name;
		stream>>name;
		name+="_";
		for (int j=0; j!=contig.size(); ++j)
		{
			stringstream stream1;
			stream1<<j;
			string name1=name+stream1.str();

//			multiCtg.appendFromFasta2(contig[j],name1);
			outs<<">"<<name1<<"\n";
			for (unsigned kk=0; kk!=contig[j].size(); ++kk)
			{
				outs<<contig[j][kk];
			}
			outs<<endl;
		}
		
		cout<<"ass() end at ends "<<i<<endl;
	
}

static void outPutAndAppCtg(std::ostream& outs){
	for (int i=0; i!=endNum;++i)
	{
		if (!CONTIG[i].empty())
		{
			stringstream stream;
			stream<<args.overlap<<'_'<<i+1;
			string name;
			stream>>name;
			name+="_";
			for (int j=0; j!=CONTIG[i].size(); ++j)
			{
				stringstream stream1;
				stream1<<j+1;
				string name1=name+stream1.str();
//				multiCtg.appendFromFasta2(CONTIG[i][j],name1);
				outs<<">"<<name1<<"\t"<<traceCount[i][j]<<"\n"<<CONTIG[i][j]<<endl;
				/*for (unsigned kk=0; kk!=CONTIG[i][j].size(); ++kk)
				{
					outs<<CONTIG[i][j][kk];
				}
				outs<<endl;*/
			}
		}
	}
}
/*
std::istream& appendFromFasta2( std::istream& in,int format, MultiSequence &multi,
				Alphabet& alph,std::size_t batchSize){
	std::size_t maxSeqBytes = batchSize;
	if( multi.finishedSequences() == 0 ) maxSeqBytes = std::size_t(-1);
	
	indexT oldUnfinishedSize = multi.unfinishedSize();
	if( format == 1 ) multi.appendFromFasta( in, maxSeqBytes );
	else if( format == 2 ) multi.appendFromFastq( in, maxSeqBytes );
	// encode the newly-read sequence
	alph.tr( multi.seqWriter() + oldUnfinishedSize,
		multi.seqWriter() + multi.unfinishedSize() );
	
	//uchar *beg=multi.seqWriter() + oldUnfinishedSize;
	//uchar *end=multi.seqWriter() + multi.unfinishedSize();
	//for(;beg<end;++beg ){
	//	std::cerr<<(int)(*beg);
	//}
	//std::cerr<<endl;

	return in;
}

void reverseComplementQuery(Alphabet& readsAlph, MultiSequence& reads){
	cerr<<"reverse complementing..."<<endl;
	readsAlph.rc( reads.seqWriter(), reads.seqWriter() + reads.finishedSize() );
//	if( format == 2 ){
//		std::reverse( reads.qualityWriter(),
//		reads.qualityWriter() +
//		reads.finishedSize() * reads.qualsPerLetter() );

	// I think the multiplication can overflow, but it's very unlikely
//	}
}
*/
void readConfig(istream& stream,vector<vector<string> >& fasta_pe,vector<vector<string> >& fastq_pe,
		vector<string>& fasta_se,vector<string>&  fastq_se,vector<unsigned>& fasta_insert,
		vector<unsigned>& fastq_insert){
	string sentence,line,word;
	vector<string> fa,fq;
	unsigned insertSize;
	string temp;
	while (getline(stream,sentence,'>')){
		if(sentence=="") continue;
		istringstream iss(sentence); //need test
		while(getline(iss,line,'\n')){
//			cout<<line<<endl;
			istringstream ist(line);
			getline(ist,word,'=');
			if(word=="insertSize") {ist>>insertSize;}
//			if(word[0]=='#') continue;  //need test
			if (word=="f1" || word=="f2"){
				ist>>temp;
				fa.push_back(temp);
				if(word=="f1") fasta_insert.push_back(insertSize);
			}
			if (word=="q1" || word=="q2"){
				ist>>temp;
//				cout<<"temp="<<temp<<endl;
				fq.push_back(temp);
				if (word=="q1") fastq_insert.push_back(insertSize);
			}
			if(word=="f"){
				ist>>temp;
				fasta_se.push_back(temp);
			}
			if(word=="q"){
				ist>>temp;
				fastq_se.push_back(temp);
			}
		}
//		cout<<"readConfig 1"<<endl;
		if(fa.size()>0) fasta_pe.push_back(fa);
//		cout<<"readConfig 2"<<endl;
		if(fq.size()>0) fastq_pe.push_back(fq);
//		cout<<"readConfig 3"<<endl;
		fa.clear();
//		cout<<"readConfig 4"<<endl;
		fq.clear();
//		cout<<"readConfig 5"<<endl;
	}
//	cout<<"readConfig 6"<<endl;
/*	std::cout<<"### 1 ###"<<"\n";
    for(indexT i=0;i<fasta_pe.size();i++){
		for(indexT j=0;j<fasta_pe[i].size();j+=2){
			std::cout<<i<<"\t"<<j<<"\t"<<fasta_pe[i][j]<<endl;
			std::cout<<i<<"\t"<<j<<"\t"<<fasta_pe[i][j+1]<<endl;
		}
	}
	for(indexT i=0;i<fastq_insert.size();i++)
	{
		std::cout<<i<<"\t"<<fastq_insert[i]<<endl;
	}
        std::cout<<"### 2 ###"<<"\n";
*/

	assert(fasta_pe.size()==fasta_insert.size());
	assert(fastq_pe.size()==fastq_insert.size());
}

void inputReads(vector<vector<string> >& fasta_pe,vector<vector<string> >& fastq_pe,
		vector<string>& fasta_se,vector<string>&  fastq_se,char flag){
	if (!fasta_pe.empty())
	{
		for (unsigned j=0; j!=fasta_pe.size(); ++j)
		{
			for (unsigned k=0; k<fasta_pe[j].size(); k+=2)
			{
				ifstream ap1,ap2;
				ap1.open(fasta_pe[j][k].c_str());
				ap2.open(fasta_pe[j][k+1].c_str());
				if (!ap1 || !ap2)
				{
					cout<<"can't open file fasta_pe[j][k] or fasta_pe[j][k+1]\n";
					exit(1);
				}
				cout<<"inputReads "<<fasta_pe[j][k]<<endl;
				cout<<"inputReads "<<fasta_pe[j][k+1]<<endl;
				readFasta(ap1,flag);  //reads number index in the layer of different files has't resove
				readFasta(ap2,flag);
				ap1.close();
				ap2.close();
				cout<<"inputReads 2 "<<endl;
			}
		}
	}
	if (!fastq_pe.empty())
	{
		for (unsigned j=0; j!=fastq_pe.size(); ++j)
		{
			for (unsigned k=0; k<fastq_pe[j].size(); k+=2)
			{
				ifstream qp1,qp2;
				qp1.open(fastq_pe[j][k].c_str());
				qp2.open(fastq_pe[j][k+1].c_str());
				readFastq(qp1,flag);
				readFastq(qp2,flag);
				qp1.close();
				qp2.close();
			}
		}
	}
	if (!fasta_se.empty())
	{
		for(unsigned j=0;j<fasta_se.size();++j){
			ifstream as;
			as.open(fasta_se[j].c_str());
			cout<<"input reads "<<fasta_se[j]<<endl;
			readFasta(as,flag);	
			as.close();
		}
	}
	if (!fastq_se.empty())
	{
		for(unsigned j=0;j<fastq_se.size();++j){
			ifstream qs;
			qs.open(fastq_se[j].c_str());
			readFastq(qs,flag);
			qs.close();
		}
	}
}

static void check1end(int t)
{
	if(flagBuf[t])
		return;
	bool late=0;
	pthread_mutex_lock(&mutex);
	if(!flagBuf[t]){
		flagBuf[t] = 1;
	}else 
		late = 1;
	pthread_mutex_unlock(&mutex);
	if(late)
		return;
//	counters[thrdID]++;
//	fill1scaf(scafCounter+t+1,ctgStackBuffer[t],thrdID);
	ass1by1(end1,end2,t);
}

static void checkWorkSignal()
{
	int t;

	while(1){
		usleep(10);
		for(t=0;t<thrd_num;t++)
			if(thrdSignal[t+1])
				break;
		if(t==thrd_num)
			break;
	}
}

static void thread_wait(pthread_t *threads)
{
	int i;
	for(i=0;i<thrd_num;i++)
		if(threads[i]!=0)
			pthread_join(threads[i],NULL);
	
}

static void *threadRoutine(void *thrID)
{
	int i;
	uchar *a;
	a=(uchar *)thrID;
	//printf("%dth thread with threadID %d, hash_table %p\n",id,prm.threadID,prm.hash_table);
	for(i=0;i<endNum;i++){
		check1end(i);
	}
	
	*a=0;
	//thrdSignal[(int)thrID+1] = 0;
}

static void creatThrds(pthread_t *threads)
{
	unsigned char i;
	int temp;
	for(i=0;i<thrd_num;i++){
		if((temp = pthread_create(&threads[i], NULL, threadRoutine, &thrdSignal[i+1])) != 0){ 
			printf("create threads failed\n");
			exit(1);
		}
	}
	printf("%d thread created\n...\n",thrd_num);
}

int main(int argc, char *argv[])
{
//	printf("Hello, world\n");
//	cout<<"1"<<endl;
//	cout<<"sizeof(ubyte8) "<<sizeof(ubyte8)<<endl;
	clock_t startTime=std::clock();
	args.fromArg(argc,argv);
//	cout<<"SOAPbarcode start time is "<<startTime/CLOCKS_PER_SEC<<"sec"<<endl;
	
	thrd_num=args.thrd_n;

	vector<vector<string> > fasta_pe,fastq_pe;
	vector<string> fasta_se,fastq_se;
	vector<unsigned> fasta_insert,fastq_insert;
	std::ifstream readStream;
	std::istream& readFilesStream = openIn( args.readsList, readStream );
//	cout<<"3"<<endl;
	readConfig(readFilesStream,fasta_pe,fastq_pe,fasta_se,fastq_se,fasta_insert,fastq_insert);
	
	ifstream ie;
	ie.open(args.endsFile.c_str());
	if (!ie)
	{
		cout<<"can't open file "<<args.endsFile<<endl;
		exit(-1);
	}

	inputEndsFile(ie,end1,end2);
	clock_t currenTime=std::clock();
        cout<<"Time spend on input file is "<<(currenTime-startTime+0.0)/CLOCKS_PER_SEC<<" sec\tcurrenTime= "<<currenTime/CLOCKS_PER_SEC<<"sec"<<endl;
        startTime=currenTime;

	//kmerul,kmerInterval,kmerInit;
	int itc=0;
	std::ofstream outStream;
	std::ostream& oStream=openOut(args.outPrefix+".contig", outStream);
	while (args.overlap>args.kmerll)
	{
		WORDFILTER=createFilter(args.overlap);
		kset.initialize(1, 0.75);//1024*1024
		inputReads(fasta_pe,fastq_pe,fasta_se,fastq_se,1);
		currenTime=std::clock();
		cout<<"Overlap="<<args.overlap<<",time spend on make graph is "<<(currenTime-startTime+0.0)/CLOCKS_PER_SEC<<" sec\tcurrenTime= "<<currenTime/CLOCKS_PER_SEC<<"sec"<<endl;
		startTime=currenTime;
		
		ubyte4 readTkmerSize=startPosition();
		readTkmer.resize(readTkmerSize,0);
		inputReads(fasta_pe,fastq_pe,fasta_se,fastq_se,0);
		resetStart();
		currenTime=std::clock();
		cout<<"Time spend on link read2kmer is "<<(currenTime-startTime+0.0)/CLOCKS_PER_SEC<<" sec\tcurrenTime= "<<currenTime/CLOCKS_PER_SEC<<"sec"<<endl;
		startTime=currenTime;

//		for (int i=0; i!=flagResult.size(); ++i)
//		{
//			cout<<i<<"\t"<<(int)flagResult[i]<<"\t"<<(int)flagBuf[i]<<endl;
//		}
//		if (itc)
//		{
//			printHash();
//			exit(1);
//		}
		if (!itc)
		{
			endNum=end1.size();
			CONTIG.resize(endNum);
			traceCount.resize(endNum);
			flagBuf.resize(endNum,0);
			flagResult.resize(endNum,0);
			itc++;
		}else{
			for (int i=0; i!=flagBuf.size(); ++i)
			{
				if (!flagResult[i])
				{
					flagBuf[i]=0;
				}
			}
		}
//		cout<<"kmer = "<<args.overlap<<endl;
//		for (int i=0; i!=flagBuf.size(); ++i)
//		{
//			cout<<i<<"\t"<<(int)flagBuf[i]<<endl;
//		}
		thrdSignal.resize(thrd_num+1,0);
		for (int i=0; i!=thrd_num; ++i)
		{
			thrdSignal[i+1]=i+1;
		}
		pthread_t threads[thrd_num];
		creatThrds(threads);
		checkWorkSignal();
		thread_wait(threads);

		kset.freeHash();
		readTkmer.clear();
		
		outPutAndAppCtg(oStream);
		CONTIG.clear();
		CONTIG.resize(end1.size());
			
		cout<<"Assembly of overlap="<<args.overlap<<" has been finished.\n";
		currenTime=std::clock();
		cout<<"Time spend on assembly is "<<(currenTime-startTime+0.0)/CLOCKS_PER_SEC<<" sec\tcurrenTime= "<<currenTime/CLOCKS_PER_SEC<<"sec"<<endl;
		startTime=currenTime;
		
		args.overlap-=args.kmerInterval;
	}
	
//	printHash();
//	exit(1);
	//add
	
//	printR2K(readTkmer);
//	exit(1);

//	removeTips();
//	currenTime=std::clock();
//	cout<<"time spend on removeTips is "<<(currenTime-startTime+0.0)/CLOCKS_PER_SEC<<" sec\tcurrenTime= "<<currenTime/CLOCKS_PER_SEC<<endl;
//    startTime=currenTime;
	
//	printHash();
//	exit(1);
//	for (unsigned i=0; i!=end1.size(); ++i)
//	{
//		cout<<end1[i]<<"\t"<<end2[i]<<endl;
//	}
//	cout<<"7"<<endl;
	
	
//	ass(end1, end2);
//	cout<<"8"<<endl;
//	printHash();
//	outPutAndAppCtg();
	CONTIG.clear();
	traceCount.clear();
	//test multiCtg;
	cout<<"The assembly has been done!"<<endl;
	return 0;
}
