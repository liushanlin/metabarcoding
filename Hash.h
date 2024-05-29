#ifndef HASH_HH
#define HASH_HH

#include "Kmer.hh"

#ifndef K_LOAD_FACTOR
#define K_LOAD_FACTOR 0.75
#endif

//#define MAX_KMER_COV 255
#define MAX_KMER_COV 65535
//#define EDGE_BIT_SIZE 6
//#define EDGE_XOR_MASK 0x3FU
#define LINKS_BITS 0x00FFFFFFU

#define get_kmer_seq(mer) ((mer).seq)
#define set_kmer_seq(mer, val) ((mer).seq = val)

#define get_kmer_num(mer) ((mer).num)
#define set_kmer_num(mer,val) ((mer).num = val)

#define get_kmer_start(mer) ((mer).start)
#define set_kmer_start(mer,val) ((mer).start = val)

#define get_kmer_left_cov1(mer, idx) (((mer).l_links1>>((idx)*EDGE_BIT_SIZE))&EDGE_XOR_MASK)
#define get_kmer_left_cov2(mer, idx) (((mer).l_links2>>((idx)*EDGE_BIT_SIZE))&EDGE_XOR_MASK)
#define set_kmer_left_cov1(mer, idx, val) ((mer).l_links1 = ((mer).l_links1&(~(EDGE_XOR_MASK<<((idx)*EDGE_BIT_SIZE)))) | (((val)&EDGE_XOR_MASK)<<((idx)*EDGE_BIT_SIZE)) )
#define set_kmer_left_cov2(mer, idx, val) ((mer).l_links2 = ((mer).l_links2&(~(EDGE_XOR_MASK<<((idx)*EDGE_BIT_SIZE)))) | (((val)&EDGE_XOR_MASK)<<((idx)*EDGE_BIT_SIZE)) )
#define get_kmer_left_covs(mer) (get_kmer_left_cov1(mer, 0) + get_kmer_left_cov1(mer, 1) + get_kmer_left_cov1(mer, 2) + get_kmer_left_cov1(mer, 3)) + (get_kmer_left_cov2(mer, 0) + get_kmer_left_cov2(mer, 1) + get_kmer_left_cov2(mer, 2) + get_kmer_left_cov2(mer, 3))

#define get_kmer_right_cov1(mer, idx) (((mer).r_links1>>((idx)*EDGE_BIT_SIZE))&EDGE_XOR_MASK)
#define get_kmer_right_cov2(mer, idx) (((mer).r_links2>>((idx)*EDGE_BIT_SIZE))&EDGE_XOR_MASK)
#define set_kmer_right_cov1(mer, idx, val) ((mer).r_links1 = ((mer).r_links1&(~(EDGE_XOR_MASK<<((idx)*EDGE_BIT_SIZE)))) | (((val)&EDGE_XOR_MASK)<<((idx)*EDGE_BIT_SIZE)) )
#define set_kmer_right_cov2(mer, idx, val) ((mer).r_links2 = ((mer).r_links2&(~(EDGE_XOR_MASK<<((idx)*EDGE_BIT_SIZE)))) | (((val)&EDGE_XOR_MASK)<<((idx)*EDGE_BIT_SIZE)) )
#define get_kmer_right_covs(mer) (get_kmer_right_cov1(mer, 0) + get_kmer_right_cov1(mer, 1) + get_kmer_right_cov1(mer, 2) + get_kmer_right_cov1(mer, 3)) + (get_kmer_right_cov2(mer, 0) + get_kmer_right_cov2(mer, 1) + get_kmer_right_cov2(mer, 2) + get_kmer_right_cov2(mer, 3))


#define is_kmer_entity_null(flags, idx)    ((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x01)
#define is_kmer_entity_del(flags, idx)     ((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x02)
#define set_kmer_entity_null(flags, idx)   ((flags)[(idx)>>4] |= (0x01u<<(((idx)&0x0f)<<1)))
#define set_kmer_entity_del(flags, idx)    ((flags)[(idx)>>4] |= (0x02u<<(((idx)&0x0f)<<1)))
#define clear_kmer_entity_null(flags, idx) ((flags)[(idx)>>4] &= ~(0x01u<<(((idx)&0x0f)<<1)))
#define clear_kmer_entity_del(flags, idx)  ((flags)[(idx)>>4] &= ~(0x02u<<(((idx)&0x0f)<<1)))
#define exists_kmer_entity(flags, idx)     (!((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x03))

namespace Barcode{
class Hash
{
	public:
		Kmer_t *array;
		ubyte4 *flags;
		ubyte8 size;
		ubyte8 count;
		ubyte8 max;
		double load_factor;
		ubyte8 iter_ptr;	

		//KmerSet* init_kmerset(ubyte8 init_size, float load_factor);
		bool KmerEqual(Kmer kmer1,Kmer kmer2);
		void initialize(ubyte8 init_size, float load_fact);
		int search_kmerset(Kmer seq, Kmer_t **rs);
		int put_kmerset(Kmer seq, ubyte left, ubyte right,Kmer_t **kmer_p);
		byte8 count_kmerset();
		void freeHash();
		void dislink2nextUncertain(Kmer_t *node,char ch,bool smaller);
		void dislink2prevUncertain(Kmer_t *node,char ch,bool smaller);
		int count_branch2prev(Kmer_t *node);
		int count_branch2next(Kmer_t *node);
		char firstCharInKmer(Kmer kmer);
//		__uint128_t Kmer2int128(Kmer seq);
		ubyte8 modular (Kmer seq);

	private:
		void update_kmer(Kmer_t *mer, ubyte left, ubyte right);
		void set_new_kmer(Kmer_t *mer, Kmer seq, ubyte left, ubyte right);
		int is_prime_kh(ubyte8 num) const;
		ubyte8 find_next_prime_kh(ubyte8 num) const;
		ubyte8 get_kmerset(Kmer seq);
		int exists_kmerset(Kmer seq);
		void encap_kmerset(ubyte8 num);
		void reset_iter_kmerset();
		ubyte8 iter_kmerset(Kmer_t **rs) ;
};
}

#endif
