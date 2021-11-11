#ifndef __QALSH_TREE_H
#define __QALSH_TREE_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>
#include <cstdint>

#include "def.h"
#include "util.h"
#include "random.h"
#include "lbtree.h"

// -----------------------------------------------------------------------------
//  Query-Aware Locality-Sensitive Hashing (QALSH) is used to solve the problem 
//  of c-Approximate Nearest Neighbor (c-ANN) search.
//
//  the idea was introduced by Qiang Huang, Jianlin Feng, Yikai Zhang, Qiong 
//  Fang, and Wilfred Ng in their paper "Query-aware locality-sensitive hashing 
//  for approximate nearest neighbor search", in Proceedings of the VLDB 
//  Endowment (PVLDB), 9(1), pages 1–12, 2015.
//
//  This is an NVM-version of QALSH, using B+-Trees in NVM
// -----------------------------------------------------------------------------
class QALSH_TREE {
public:
	void check();
	// -------------------------------------------------------------------------
	// QALSH_TREE(							// constructor
	// 	int   n,						// cardinality
	// 	int   d,						// dimensionality
	// 	float p,						// l_p distance
	// 	float zeta,						// symmetric factor of p-stable distr.
	// 	float ratio);					// approximation ratio

	// -------------------------------------------------------------------------
	QALSH_TREE(							// constructor
		int   n,						// cardinality
		int   d,						// dimensionality
		float p,						// l_p distance
		float zeta,						// symmetric factor of p-stable distr.
		float ratio,					// approximation ratio
		const char* nvm_file,           // nvm file
		const char* data_file,          // data file on disk
		const float** data);			// data objects
	
	// -------------------------------------------------------------------------
	/**
	 * data on disk
	*/
	// QALSH_TREE(							// constructor
	// 	int   n,						// cardinality
	// 	int   d,						// dimensionality
	// 	float p,						// l_p distance
	// 	float zeta,						// symmetric factor of p-stable distr.
	// 	float ratio,					// approximation ratio
	// 	int group_size,
	// 	const char* nvm_file,           // nvm file
	// 	const char* data_file);		    // data file on disk
	
	// -------------------------------------------------------------------------
	/**
	 * Restore from NVM
	*/
	QALSH_TREE(
		const char* data_file_,
		const char* nvm_file
	);

	// -------------------------------------------------------------------------
	~QALSH_TREE();						// destructor

	// -------------------------------------------------------------------------
	float calc_hash_value(			// calc hash value
		int   tid,						// hash table id
		const float *data);				// one data/query object

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	// int knn(						// k-NN search
	// 	int   top_k,					// top-k value
	// 	const float *query,				// input query object
	// 	MinK_List *list);				// k-NN results (return)

	// -------------------------------------------------------------------------
	int knn(						// k-NN search
		int   top_k,					// top-k value
		const float *query,				// input query object
		std::vector<KV> &cand);		// k-NN results (return)

	// -------------------------------------------------------------------------
	// int knn(						// k-NN search
	// 	int   top_k,					// top-k value
	// 	float R,						// limited search range
	// 	const float *query,				// input query object
	// 	std::vector<int> &cand);		// object id mapping

	void getDataById(int id, float* data);

	// -------------------------------------------------------------------------
	int    n_pts_;					// cardinality
	int    dim_;					// dimensionality
	float  p_;						// l_p distance
	float  zeta_;					// symmetric factor of p-stable distr.
	float  ratio_;					// approximation ratio
	float  w_;						// bucket width
	int    m_;						// number of hashtables
	int    l_;						// collision threshold

	float          **a_;			// hash functions
	tree*          *tables_;		// hash tables
	FILE*          fp_;             // data file pointer

	// NVM related
	// data in NVM: 256 header(Root)
	char*          nvm_pool_;

protected:
	// -------------------------------------------------------------------------
	void init();					// basic initialzation

	void bulkload(const float** data);

	// -------------------------------------------------------------------------
	float calc_l0_prob(				// calc <p1> and <p2> for L_{0.5} distance
		float x);						// x = w / (2.0 * r)

	float calc_l1_prob(				// calc <p1> and <p2> for L_{1.0} distance
		float x);						// x = w / (2.0 * r)

	float calc_l2_prob(				// calc <p1> and <p2> for L_{2.0} distance
		float x);						// x = w / (2.0 * r)
};

#endif // __QALSH_TREE_H
