#ifndef KMER_HH
#define KMER_HH

#include "Type.hh"

//#define EDGE_BIT_SIZE 6
//#define EDGE_BIT_SIZE 7
//#define EDGE_BIT_SIZE 8
#define EDGE_BIT_SIZE 16
//#define EDGE_XOR_MASK 0x3FU
//#define EDGE_XOR_MASK 0x7FU
//#define EDGE_XOR_MASK 0xFFU
#define EDGE_XOR_MASK 0xFFFFU
#define base2int(base) (char)(((base)&0x06)>>1)
#define int2base(seq) "ACTG"[seq]
#define int2compbase(seq) "TGAC"[seq]
#define int_comp(seq) (char)(seq^0x02) //(char)((0x4E>>((seq)<<1))&0x03)

namespace Barcode{
/*typedef struct Kmer
{
	unsigned long long high,low;
};
*/
typedef struct Kmer
{
	unsigned long long high1, low1, high2, low2;
};

typedef struct Kmer_t
{
	Kmer seq;
	ubyte4 l_links1;
	ubyte4 l_links2;
	ubyte4 r_links1;
	ubyte4 r_links2;
	ubyte4 start;
	ubyte4 num;
};

/*typedef struct Kmer_t
{
	Kmer seq;
//	ubyte4 l_links;
	ubyte4 l_links:4*EDGE_BIT_SIZE;
	ubyte4 r_links:4*EDGE_BIT_SIZE;
	ubyte4 linear:1;
	ubyte4 deleted:1;
	ubyte4 checked:1;
	ubyte4 single:1;
	ubyte4 twin:2;
	ubyte4 inEdge:2;
};*/
}

#endif
