#include <math.h>
#include <iostream>
#include "Hash.hh"

using namespace Barcode;

//static const Kmer_t empty_kmer = {{0,0}, 0, 0, 0, 0, 0, 1, 0,0};
//static const Kmer_t empty_kmer = {{0,0,0,0}, 0, 0, 0, 0, 0, 1, 0,0};
static const Kmer_t empty_kmer = {{0,0,0,0}, 0, 0,0,0};
/*
__uint128_t Hash::Kmer2int128(Kmer seq)
{
	__uint128_t temp;
	temp = seq.high;
	temp <<= 64;
	temp |= seq.low;

	return temp;
}*/

ubyte8 Hash::modular (Kmer seq)
{
	ubyte8 temp;

	temp = (seq.high1 % size) << 32 | (seq.low1 >> 32 & 0xffffffff);
	temp = (temp % size) << 32 | (seq.low1 & 0xffffffff);
	temp = (temp % size) << 32 | (seq.high2 >> 32 & 0xffffffff);
	temp = (temp % size) << 32 | (seq.high2 & 0xffffffff);
	temp = (temp % size) << 32 | (seq.low2 >> 32 & 0xffffffff);
	temp = (temp % size) << 32 | (seq.low2 & 0xffffffff);
	temp = (ubyte8) (temp % size);
	return temp;
}
/*
bool Hash::KmerEqual(Kmer kmer1,Kmer kmer2)
{
	if(kmer1.high==kmer2.high&&kmer1.low==kmer2.low)
		return 1;
	else return 0;
}*/
bool Hash::KmerEqual (Kmer kmer1, Kmer kmer2)
{
	if (kmer1.low2 != kmer2.low2 || kmer1.high2 != kmer2.high2 || kmer1.low1 != kmer2.low1 || kmer1.high1 != kmer2.high1)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

void Hash::update_kmer(Kmer_t *mer, ubyte left, ubyte right){
	ubyte4 cov;

	if(left<4){
		if (left<2)
		{
			cov = get_kmer_left_cov1(*mer, left);
			if(cov < MAX_KMER_COV){
				set_kmer_left_cov1(*mer, left, cov + 1);
			}	
		}else{
			cov = get_kmer_left_cov2(*mer, left-2);
			if(cov < MAX_KMER_COV){
				set_kmer_left_cov2(*mer, left-2, cov + 1);
			}
		}
//		std::cout<<"update_kmer left "<<(int)left<<"\t"<<cov<<std::endl;
//		cov = get_kmer_left_cov2(*mer, left);
//		std::cout<<"* update_kmer left "<<(int)left<<"\t"<<cov<<std::endl;
	}

	if(right<4){
		if (right<2)
		{
			cov = get_kmer_right_cov1(*mer, right);
			if(cov < MAX_KMER_COV){
				set_kmer_right_cov1(*mer, right, cov + 1);
			}	
		}else{
			cov = get_kmer_right_cov2(*mer, right-2);
			if(cov < MAX_KMER_COV){
				set_kmer_right_cov2(*mer, right-2, cov + 1);
			}
		}

//		std::cout<<"update_kmer right "<<(int)right<<"\t"<<cov<<std::endl;
//		cov = get_kmer_right_cov1(*mer, right);
//		std::cout<<"* update_kmer right "<<(int)right<<"\t"<<cov<<std::endl;

	}
}

void Hash::set_new_kmer(Kmer_t *mer, Kmer seq, ubyte left, ubyte right){
	*mer = empty_kmer;
	set_kmer_seq(*mer, seq);
	if(left<4){
		if (left<2)
		{
			set_kmer_left_cov1(*mer, left, 1);
		}else{
			set_kmer_left_cov2(*mer, left-2, 1);
		}
	}

//	std::cout<<"left "<<(int)left<<"\t"<<1<<std::endl;
//	ubyte4 cov = get_kmer_left_cov2(*mer, left);
//	std::cout<<"* left "<<(int)left<<"\t"<<cov<<std::endl;
	
	if(right<4){
		if (right<2)
		{
			set_kmer_right_cov1(*mer, right, 1);
		}else{
			set_kmer_right_cov2(*mer, right-2, 1);
		}
	}

//	std::cout<<"right "<<(int)right<<"\t"<<1<<std::endl;
//	cov = get_kmer_right_cov1(*mer, right);
//	std::cout<<"* right "<<(int)right<<"\t"<<cov<<std::endl;

}


int Hash::is_prime_kh(ubyte8 num) const{
	ubyte8 i, max;
	if(num < 4) return 1;
	if(num % 2 == 0) return 0;
	max = (ubyte8)sqrt((float)num);
	for(i=3;i<max;i+=2){ if(num % i == 0) return 0; }
	return 1;
}

ubyte8 Hash::find_next_prime_kh(ubyte8 num) const{
	if(num % 2 == 0) num ++;
	while(1){ if(is_prime_kh(num)) return num; num += 2; }
}

void Hash::initialize(ubyte8 init_size, float load_fact){
//	KmerSet *set;
	if(init_size < 3) init_size = 3;
	else init_size = find_next_prime_kh(init_size);
//	set = (KmerSet*)malloc(sizeof(KmerSet));
	size   = init_size;
	count  = 0;
	max    = ubyte8(size * load_fact);
	if(load_fact <= 0) load_fact = 0.25f;
	else if(load_fact >= 1) load_fact = 0.75f;
	load_factor = load_fact;
	iter_ptr    = 0;
	array = (Kmer_t *)calloc(size, sizeof(Kmer_t));
	flags = (ubyte4 *)malloc((size + 15)/16 * 4);
	memset(flags, 0x55, (size + 15) / 16 * 4);
//	return set;
}

ubyte8 Hash::get_kmerset(Kmer seq) {
	ubyte8 hc;
//	__uint128_t temp;
//	temp = Kmer2int128(seq);
//	hc = temp % size;
	hc = modular (seq);
	while(1){
		if(is_kmer_entity_null(flags, hc)){
			return hc;
		} else {
			if(KmerEqual(get_kmer_seq(array[hc]),seq)) return hc;
		}
		hc ++;
		if(hc == size) hc = 0;
	}
	return size;
}

int Hash::search_kmerset(Kmer seq, Kmer_t **rs){
	ubyte8 hc;
//	__uint128_t temp;
//	temp = Kmer2int128(seq);
//	hc = temp % size;
	hc = modular (seq);
	while(1){
		if(is_kmer_entity_null(flags, hc)){
			return 0;
		} else {
			if(KmerEqual(get_kmer_seq(array[hc]),seq)){
				*rs = array + hc;
				return 1;
			}
		}
		hc ++;
		if(hc == size) hc = 0;
	}
	return 0;
}

int Hash::exists_kmerset(Kmer seq){ 
	ubyte8 idx;
	idx = get_kmerset(seq);
	return !is_kmer_entity_null(flags, idx);
}

void Hash::encap_kmerset(ubyte8 num){
	ubyte4 *flags1, *f;
	ubyte8 i, n, size1, hc;
	Kmer_t key, tmp;
	if(count + num <= max) return;
//	std::cout<<"encap_kmerset 1\t"<<size<<"\t"<<max<<"\t"<<count<<std::endl;
	n = size;
	do{
		if(n < 0xFFFFFFFU)
			n <<= 1;
		else
			n += 0xFFFFFFU;
		n = find_next_prime_kh(n); 
//		std::cout<<n<<"\t"<<load_factor<<"\t"<<count<<"\t"<<num<<std::endl;
	} while(n * load_factor < count + num);
//	std::cout<<"encap_kmerset 2\t"<<size<<"\t"<<max<<"\t"<<count<<std::endl;
	array = (Kmer_t *)realloc(array, n * sizeof(Kmer_t));
	//printf("Allocate Mem %lld(%d*%lld*%d)bytes\n",thrd_num*n*sizeof(Kmer_t),thrd_num,n,sizeof(Kmer_t));
	if(array == NULL){
		fprintf(stderr, "-- Out of memory --\n");
		abort();
	}
//	std::cout<<"encap_kmerset 3\t"<<size<<"\t"<<max<<"\t"<<count<<std::endl;
	flags1 = (ubyte4 *)malloc((n+15)/16 * 4);
	memset(flags1, 0x55, (n+15)/16 * 4);
	size1 = size;
	size = n;
	max =ubyte8(n * load_factor);
	f = flags;
	flags = flags1;
	flags1 = f;
//	__uint128_t temp;
//	std::cout<<"encap_kmerset 4\t"<<size<<"\t"<<max<<"\t"<<count<<std::endl;
	for(i=0;i<size1;i++){
		if(!exists_kmer_entity(flags1, i)) continue;
		key = array[i];
		set_kmer_entity_del(flags1, i);
		while(1){
//			temp = Kmer2int128(get_kmer_seq(key));
//			hc = temp % size;
			hc = modular (get_kmer_seq (key));
			while(!is_kmer_entity_null(flags, hc)){ hc ++; if(hc == size) hc = 0; }
			clear_kmer_entity_null(flags, hc);
			if(hc < size1 && exists_kmer_entity(flags1, hc)){ 
				tmp = key;
				key = array[hc];
				array[hc] = tmp;
				set_kmer_entity_del(flags1, hc);
			} else {
				array[hc] = key;
				break;
			}
		}
	}
	free(flags1);
//	std::cout<<"encap_kmerset 5\t"<<size<<"\t"<<max<<"\t"<<count<<std::endl;
}

int Hash::put_kmerset(Kmer seq, ubyte left, ubyte right, Kmer_t **kmer_p){
	ubyte8 hc;
	encap_kmerset(1);
//	std::cout<<"put_kmerset "<<size<<"\t"<<max<<"\t"<<count<<std::endl;
//	__uint128_t temp;
//	temp = Kmer2int128(seq);
//	hc = temp % size;
	hc = modular (seq);
	do{
		if(is_kmer_entity_null(flags, hc)){
			clear_kmer_entity_null(flags, hc);
			set_new_kmer(array + hc, seq, left, right);
			set_kmer_num(array[hc],1);//+
			count ++;
			*kmer_p = array + hc;
			return 0;
		} else {
			if(KmerEqual(get_kmer_seq(array[hc]),seq)){
				update_kmer(array + hc, left, right);
				ubyte4 cov=get_kmer_num(array[hc]); //+
				set_kmer_num(array[hc],cov+1);//+
//				array[hc].single = 0;
				*kmer_p = array + hc;
				return 1;
			}
		}
		hc ++;
		if(hc == size) hc = 0;
	} while(1);
	*kmer_p = NULL;
	return 0;
}

byte8 Hash::count_kmerset(){ return count; }

void Hash::reset_iter_kmerset(){ iter_ptr = 0; }

ubyte8 Hash::iter_kmerset(Kmer_t **rs){
	while(iter_ptr < size){
		if(!is_kmer_entity_null(flags, iter_ptr)){
			*rs = array + iter_ptr;
			iter_ptr ++;
			return 1;
		}
		iter_ptr ++;
	}
	return 0;
}

void Hash::freeHash()
{
        free((Kmer_t *)array);
        free((ubyte4 *)flags);
}

/*
void Hash::free_kmerset(KmerSet *set){
	free(array);
	free(flags);
	free(set);
}

void Hash::free_Sets(KmerSet **sets,int num){
	int i;
	for(i=0;i<num;i++)
		free_kmerset(sets[i]);
	free((void*)sets);
}*/

int Hash::count_branch2prev(Kmer_t *node)
{
	int num = 0,i;

	for(i=0;i<4;i++){
		if (i<2)
		{
			if(get_kmer_left_cov1(*node,i)>0)
				num++;
		}else{
			if(get_kmer_left_cov2(*node,i-2)>0)
				num++;
		}
		
	}
	return num;
}

int Hash::count_branch2next(Kmer_t *node)
{
	int num = 0,i;

	for(i=0;i<4;i++){
		if (i<2)
		{
			if(get_kmer_right_cov1(*node,i)>0)
				num++;
		}else{
			if(get_kmer_right_cov2(*node,i-2)>0)
				num++;
		}
	}
	return num;
}

void Hash::dislink2prevUncertain(Kmer_t *node,char ch,bool smaller)
{
	if(smaller){
		if (ch<2)
		{
			set_kmer_left_cov1(*node,ch,0);	
		}else{
			set_kmer_left_cov2(*node,ch-2,0);
		}
	}else{
		if (int_comp(ch)<2)
		{
			set_kmer_right_cov1(*node,int_comp(ch),0);
		}else{
			set_kmer_right_cov2(*node,int_comp(ch)-2,0);
		}
	}
}

void Hash::dislink2nextUncertain(Kmer_t *node,char ch,bool smaller)
{
	if(smaller){
		if (ch<2)
		{
			set_kmer_right_cov1(*node,ch,0);
		}else{
			set_kmer_right_cov2(*node,ch-2,0);
		}
	}else{
		if (int_comp(ch)<2)
		{
			set_kmer_left_cov1(*node,int_comp(ch),0);
		}else{
			set_kmer_left_cov2(*node,int_comp(ch)-2,0);
		}
	}
}
	
