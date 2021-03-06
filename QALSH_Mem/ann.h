#ifndef __ANN_H
#define __ANN_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <sys/time.h>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "qalsh.h"
#include "qalsh_tree.h"

struct Result;

// -----------------------------------------------------------------------------
int linear_scan(					// k-NN search by linear scan
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int qalsh(							// k-NN search by qalsh
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	int   k,
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *out_path);				// output path
	
// -----------------------------------------------------------------------------
int qalsh_tree(							// k-NN search by qalsh
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	int   k,
	const float **data,					// data set
	const char*   data_file,            // data_file
	const char*   nvm_file,             // nvm file
	const float **query,				// query set
	const Result **R,					// truth set
	const char *out_path);				// output path

#endif // __ANN_H
