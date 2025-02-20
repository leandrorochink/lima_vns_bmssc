#ifndef RANDOM_H_
#define RANDOM_H_

#include <random>

using namespace std;

class Random {

private:
	double seed;

public:
	Random(int _seed);
	double getSeed();

	int get_rand_ij(int i, int j );
	int get_rand(int size );
	double get_rand01();
	double randp();
	int trand();

	template<typename _RandomAccessIterator>
	inline void random_shuffle(_RandomAccessIterator __first, _RandomAccessIterator __last){

	    if (__first != __last){
	        for (_RandomAccessIterator __i = __first + 1; __i != __last; ++__i){
	            _RandomAccessIterator __j = __first
	                    + get_rand_ij(0 , ((__i - __first)));
	            if (__i != __j)
	                std::iter_swap(__i, __j);
	        }
	    }
	}

};

#endif /* RANDOM_H_ */
