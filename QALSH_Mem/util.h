#ifndef __UTIL_H
#define __UTIL_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>
#include <chrono>
#include <ratio>

#include <unistd.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>

#include "def.h"
// #include "util.h"
#include "pri_queue.h"

struct Result;
class  MinK_List;

extern timeval g_start_time;		// global parameter: start time
extern timeval g_end_time;			// global parameter: end time

extern uint64_t g_memory;			// global parameter: memory usage
extern uint64_t g_memory_nvm;		// global parameter: memory usage in NVM
extern float    g_runtime;			// global parameter: running time
extern float    g_ratio;			// global parameter: overall ratio
extern float    g_recall;			// global parameter: recall
extern std::chrono::high_resolution_clock::time_point g_start_time2;
extern std::chrono::high_resolution_clock::time_point g_end_time2;
extern double   g_runtime2;


struct KV {
	float    key_;              // 4B key
	uint32_t id_;	     		// 4B id
		bool operator<(const KV& other) const {
		bool ret = 0;
		if (this->key_ < other.key_) {
			ret = 1;
		} 
		else if (this->key_ > other.key_) {
			ret = 0;
		} 
		else {
			if (this->id_ < other.id_) ret = 1;
			else if (this->id_ > other.id_) ret = 0;
		}
		return ret;
	}

	bool operator<(const float& other) const {
		bool ret = 0;
		if (this->key_ < other) {
			ret = 1;
		} 
		else if (this->key_ > other) {
			ret = 0;
		}
		return ret;
	}

	bool operator>(const KV& other) const {
		bool ret = 0;
		if (this->key_ < other.key_) {
			ret = 0;
		} 
		else if (this->key_ > other.key_) {
			ret = 1;
		} 
		else {
			if (this->id_ < other.id_) ret = 0;
			else if (this->id_ > other.id_) ret = 1;
		}
		return ret;
	}

	bool operator>(const float& other) const {
		bool ret = 0;
		if (this->key_ < other) {
			ret = 0;
		} 
		else if (this->key_ > other) {
			ret = 1;
		}
		return ret;
	}
};

struct Root {
	uint32_t    n_pts_;			// 4B cardinality
	uint32_t    dim_;			// 4B dimensionality
	float       p_;				// 4B l_p distance
	float       zeta_;			// 4B symmetric factor of p-stable distr.
	float       ratio_;			// 4B approximation ratio
	float       w_;				// 4B bucket width
	uint32_t    m_;				// 4B number of hashtables
	uint32_t    l_;				// 4B collision threshold
};

// -----------------------------------------------------------------------------
int KVComp(						// compare function for qsort (ascending)
	const void *e1,						// 1st element
	const void *e2);					// 2nd element

// -----------------------------------------------------------------------------
void create_dir(					// create dir if the path exists
	char *path);						// input path

// -----------------------------------------------------------------------------
int read_txt_data(					// read data (text) from disk
	int   n,							// number of data/query objects
	int   d,			 				// dimensionality
	const char *fname,					// address of data/query set
	float **data);						// data/query objects (return)

// -----------------------------------------------------------------------------
int read_bin_data(					// read data (binary) from disk
	int   n,							// number of data/query objects
	int   d,			 				// dimensionality
	const char *fname,					// address of data/query set
	float **data);						// data/query objects (return)

// -----------------------------------------------------------------------------
int read_ground_truth(				// read k-NN ground truth 
	int qn,								// number of query objects
	const char *fname,					// address of truth set
	Result **R);						// ground truth results (return)

// -----------------------------------------------------------------------------
float calc_lp_dist(					// calc L_{p} norm
	int   dim,							// dimension
	float p,							// the p value of Lp norm, p in (0, 2]
	float threshold,					// threshold
	const float *vec1,					// 1st point
	const float *vec2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_l2_sqr(					// calc l2 square distance
	int   dim,							// dimension
	float threshold,					// threshold
	const float *p1,					// 1st point
	const float *p2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_l1_dist(					// calc Manhattan distance
	int   dim,							// dimension
	float threshold,					// threshold
	const float *p1,					// 1st point
	const float *p2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_l0_sqrt(					// calc L_{0.5} sqrt distance
	int   dim,							// dimension
	float threshold,					// threshold
	const float *p1,					// 1st point
	const float *p2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_lp_pow(					// calc L_p pow_p distance
	int   dim,							// dimension
	float p,							// the p value of Lp norm, p in (0,2]
	float threshold,					// threshold
	const float *p1,					// 1st point
	const float *p2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_inner_product(			// calc inner product
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	MinK_List *list);					// results returned by algorithms

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	const Result *result);				// results returned by algorithms

float calc_recall(					// calc recall (percentage)
	int  k,								// top-k value
	const Result *R,					// ground truth results 
	std::vector<KV>& cand);		// results returned by algorithms

// -----------------------------------------------------------------------------
int ground_truth(					// find ground truth
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set);			// address of truth set

// -----------------------------------------------------------------------------
void k_nn_search(					// k-NN search
	int   n, 							// cardinality
	int   qn,							// query number
	int   d, 							// dimensionality
	int   k,							// top-k value
	float p,							// the p value of Lp norm, p in (0,2]
	const float **data,					// data objects
	const float **query,				// query objects
	Result **result);					// k-MIP results (return)

// -----------------------------------------------------------------------------
void k_nn_search(					// k-NN search
	int   n, 							// cardinality
	int   qn,							// query number
	int   d, 							// dimensionality
	int   k,							// top-k value
	float p,							// the p value of Lp norm, p in (0,2]
	const char *data,					// data and query file on NVM
	Result **result);					// k-MIP results (return)

void write_nvm_data(
	int   n,
	int   qn,
	int   d,
	const float **data,					// data objects
	const float **query,				// query objects
	const char* file_path);

void test_rebuild_tree(					// k-NN search by qalsh
	int   n,							// number of data  objects
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	const float **data,					// data set
	const char*   data_file,            // data_file
	const char*   nvm_file,             // nvm file
	const char *out_path); 				// output path

#endif // __UTIL_H
