#include "util.h"
#include "qalsh_tree.h"
#include <libpmem.h>

timeval  g_start_time;
timeval  g_end_time;

uint64_t g_memory  = 0;
uint64_t g_memory_nvm = 0;
float    g_runtime = -1.0f;
float    g_ratio   = -1.0f;
float    g_recall  = -1.0f;
std::chrono::high_resolution_clock::time_point g_start_time2;
std::chrono::high_resolution_clock::time_point g_end_time2;
double   g_runtime2 = 0.0;

// -----------------------------------------------------------------------------
void create_dir(					// create dir if the path exists
	char *path)							// input path
{
	int len = (int) strlen(path);
	for (int i = 0; i < len; ++i) {
		if (path[i] == '/') {
			char ch = path[i + 1];
			path[i + 1] = '\0';
									// check whether the directory exists
			int ret = access(path, F_OK);
			if (ret != 0) {			// create the directory
				ret = mkdir(path, 0755);
				if (ret != 0) printf("Could not create %s\n", path);
			}
			path[i + 1] = ch;
		}
	}
}

// -----------------------------------------------------------------------------
int read_txt_data(					// read data (text) from disk
	int   n,							// number of data/query objects
	int   d,			 				// dimensionality
	const char *fname,					// address of data/query set
	float **data)						// data/query objects (return)
{
	gettimeofday(&g_start_time, NULL);
	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s.\n", fname);
		return 1;
	}

	int i = 0;
	int j = 0;
	while (!feof(fp) && i < n) {
		fscanf(fp, "%d", &j);
		for (j = 0; j < d; ++j) {
			fscanf(fp, " %f", &data[i][j]);
		}
		fscanf(fp, "\n");
		++i;
	}
	assert(feof(fp) && i == n);
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float running_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Read Data: %f Seconds\n\n", running_time);

	return 0;
}

// -----------------------------------------------------------------------------
int read_bin_data(					// read data (binary) from disk
	int   n,							// number of data points
	int   d,							// dimensionality
	const char *fname,					// address of data
	float **data)						// data/query objects (return)
{
	gettimeofday(&g_start_time, NULL);
	FILE *fp = fopen(fname, "rb");
	if (!fp) {
		printf("Could not open %s\n", fname);
		return 1;
	}

	int i = 0;
	while (!feof(fp) && i < n) {
		fread(data[i++], SIZEFLOAT, d, fp);
	}
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float running_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Read Data: %f Seconds\n\n", running_time);

	return 0;
}

// -----------------------------------------------------------------------------
int read_ground_truth(				// read ground truth results from disk
	int qn,								// number of query objects
	const char *fname,					// address of truth set
	Result **R)							// ground truth results (return)
{
	gettimeofday(&g_start_time, NULL);
	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s\n", fname);
		return 1;
	}

	int tmp1 = -1;
	int tmp2 = -1;
	fscanf(fp, "%d %d\n", &tmp1, &tmp2);
	assert(tmp1 == qn && tmp2 == MAXK);

	for (int i = 0; i < qn; ++i) {
		for (int j = 0; j < MAXK; ++j) {
			fscanf(fp, "%d %f ", &R[i][j].id_, &R[i][j].key_);
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	gettimeofday(&g_end_time, NULL);
	float running_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Read Ground Truth: %f Seconds\n\n", running_time);

	return 0;
}

// -----------------------------------------------------------------------------
float calc_lp_dist(					// calc L_{p} norm
	int   dim,							// dimension
	float p,							// the p value of Lp norm, p in (0,2]
	float threshold,					// threshold
	const float *vec1,					// 1st point
	const float *vec2)					// 2nd point
{
	if (fabs(p - 2.0f) < FLOATZERO) {
		return sqrt(calc_l2_sqr(dim, SQR(threshold), vec1, vec2));
	}
	else if (fabs(p - 1.0f) < FLOATZERO) {
		return calc_l1_dist(dim, threshold, vec1, vec2);
	}
	else if (fabs(p - 0.5f) < FLOATZERO) {
		float ret = calc_l0_sqrt(dim, sqrt(threshold), vec1, vec2);
		return SQR(ret);
	}
	else {
		float ret = calc_lp_pow(dim, p, pow(threshold, p), vec1, vec2);
		return pow(ret, 1.0f / p);
	}
}

// -----------------------------------------------------------------------------
float calc_l2_sqr(					// calc l2 square distance
	int   dim,							// dimension
	float threshold,					// threshold
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	unsigned d = dim & ~unsigned(7);
	const float *aa = p1, *end_a = aa + d;
	const float *bb = p2, *end_b = bb + d;

	float r = 0.0f;
	float r0, r1, r2, r3, r4, r5, r6, r7;

	const float *a = end_a, *b = end_b;

	r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0.0f;
	switch (dim & 7) {
		case 7: r6 = SQR(a[6] - b[6]);
		case 6: r5 = SQR(a[5] - b[5]);
		case 5: r4 = SQR(a[4] - b[4]);
		case 4: r3 = SQR(a[3] - b[3]);
		case 3: r2 = SQR(a[2] - b[2]);
		case 2: r1 = SQR(a[1] - b[1]);
		case 1: r0 = SQR(a[0] - b[0]);
	}

	a = aa; b = bb;
	for (; a < end_a; a += 8, b += 8) {
		r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
		if (r > threshold) return r;

		r0 = SQR(a[0] - b[0]);
		r1 = SQR(a[1] - b[1]);
		r2 = SQR(a[2] - b[2]);
		r3 = SQR(a[3] - b[3]);
		r4 = SQR(a[4] - b[4]);
		r5 = SQR(a[5] - b[5]);
		r6 = SQR(a[6] - b[6]);
		r7 = SQR(a[7] - b[7]);
	}
	r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
	
	return r;
}

// -----------------------------------------------------------------------------
float calc_l1_dist(					// calc Manhattan distance
	int   dim,							// dimension
	float threshold,					// threshold
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	unsigned d = dim & ~unsigned(7);
	const float *aa = p1, *end_a = aa + d;
	const float *bb = p2, *end_b = bb + d;

	float r = 0.0f;
	float r0, r1, r2, r3, r4, r5, r6, r7;

	const float *a = end_a, *b = end_b;

	r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0.0f;
	switch (dim & 7) {
		case 7: r6 = fabs(a[6] - b[6]);
		case 6: r5 = fabs(a[5] - b[5]);
		case 5: r4 = fabs(a[4] - b[4]);
		case 4: r3 = fabs(a[3] - b[3]);
		case 3: r2 = fabs(a[2] - b[2]);
		case 2: r1 = fabs(a[1] - b[1]);
		case 1: r0 = fabs(a[0] - b[0]);
	}

	a = aa; b = bb;
	for (; a < end_a; a += 8, b += 8) {
		r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
		if (r > threshold) return r;

		r0 = fabs(a[0] - b[0]);
		r1 = fabs(a[1] - b[1]);
		r2 = fabs(a[2] - b[2]);
		r3 = fabs(a[3] - b[3]);
		r4 = fabs(a[4] - b[4]);
		r5 = fabs(a[5] - b[5]);
		r6 = fabs(a[6] - b[6]);
		r7 = fabs(a[7] - b[7]);
	}
	r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
	
	return r;
}

// -----------------------------------------------------------------------------
float calc_l0_sqrt(					// calc L_{0.5} sqrt distance
	int   dim,							// dimension
	float threshold,					// threshold
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	unsigned d = dim & ~unsigned(7);
	const float *aa = p1, *end_a = aa + d;
	const float *bb = p2, *end_b = bb + d;

	float r = 0.0f;
	float r0, r1, r2, r3, r4, r5, r6, r7;

	const float *a = end_a, *b = end_b;

	r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0.0f;
	switch (dim & 7) {
		case 7: r6 = sqrt(fabs(a[6] - b[6]));
		case 6: r5 = sqrt(fabs(a[5] - b[5]));
		case 5: r4 = sqrt(fabs(a[4] - b[4]));
		case 4: r3 = sqrt(fabs(a[3] - b[3]));
		case 3: r2 = sqrt(fabs(a[2] - b[2]));
		case 2: r1 = sqrt(fabs(a[1] - b[1]));
		case 1: r0 = sqrt(fabs(a[0] - b[0]));
	}

	a = aa; b = bb;
	for (; a < end_a; a += 8, b += 8) {
		r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
		if (r > threshold) return r;

		r0 = sqrt(fabs(a[0] - b[0]));
		r1 = sqrt(fabs(a[1] - b[1]));
		r2 = sqrt(fabs(a[2] - b[2]));
		r3 = sqrt(fabs(a[3] - b[3]));
		r4 = sqrt(fabs(a[4] - b[4]));
		r5 = sqrt(fabs(a[5] - b[5]));
		r6 = sqrt(fabs(a[6] - b[6]));
		r7 = sqrt(fabs(a[7] - b[7]));
	}
	r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
	
	return r;
}

// -----------------------------------------------------------------------------
float calc_lp_pow(					// calc L_p pow_p distance
	int   dim,							// dimension
	float p,							// the p value of L_p norm, p in (0,2]
	float threshold,					// threshold
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	unsigned d = dim & ~unsigned(7);
	const float *aa = p1, *end_a = aa + d;
	const float *bb = p2, *end_b = bb + d;

	float r = 0.0f;
	float r0, r1, r2, r3, r4, r5, r6, r7;

	const float *a = end_a, *b = end_b;

	r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0.0f;
	switch (dim & 7) {
		case 7: r6 = pow(fabs(a[6] - b[6]), p);
		case 6: r5 = pow(fabs(a[5] - b[5]), p);
		case 5: r4 = pow(fabs(a[4] - b[4]), p);
		case 4: r3 = pow(fabs(a[3] - b[3]), p);
		case 3: r2 = pow(fabs(a[2] - b[2]), p);
		case 2: r1 = pow(fabs(a[1] - b[1]), p);
		case 1: r0 = pow(fabs(a[0] - b[0]), p);
	}

	a = aa; b = bb;
	for (; a < end_a; a += 8, b += 8) {
		r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
		if (r > threshold) return r;

		r0 = pow(fabs(a[0] - b[0]), p);
		r1 = pow(fabs(a[1] - b[1]), p);
		r2 = pow(fabs(a[2] - b[2]), p);
		r3 = pow(fabs(a[3] - b[3]), p);
		r4 = pow(fabs(a[4] - b[4]), p);
		r5 = pow(fabs(a[5] - b[5]), p);
		r6 = pow(fabs(a[6] - b[6]), p);
		r7 = pow(fabs(a[7] - b[7]), p);
	}
	r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
	
	return r;
}

// -----------------------------------------------------------------------------
float calc_inner_product(			// calc inner product
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	float r = 0.0f;
	for (int i = 0; i < dim; ++i) {
		r += p1[i] * p2[i];
	}
	return r;
}

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int  k,								// top-k value
	const Result *R,					// ground truth results 
	MinK_List *list)					// results returned by algorithms
{
	int i = k - 1;
	int last = k - 1;
	while (i >= 0 && list->ith_key(i) > R[last].key_) {
		i--;
	}
	return (i + 1) * 100.0f / k;
}

float calc_recall(					// calc recall (percentage)
	int  k,								// top-k value
	const Result *R,					// ground truth results 
	std::vector<KV>& cand)			// results returned by algorithms
{
	int i = k - 1;
	int last = k - 1;
	while (i >= 0 && cand[i].key_ > R[last].key_) {
		i--;
	}
	return (i + 1) * 100.0f / k;
}

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int k,								// top-k value
	const Result *R,					// ground truth results 
	const Result *result)				// results returned by algorithms
{
	int i = k - 1;
	int last = k - 1;
	while (i >= 0 && result[i].key_ > R[last].key_) {
		i--;
	}
	return (i + 1) * 100.0f / k;
}

// -----------------------------------------------------------------------------
int ground_truth(					// find ground truth
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	const float **data,					// data set
	const float **query,				// query set
	const char  *truth_set)				// address of truth set
{
	gettimeofday(&g_start_time, NULL);
	FILE *fp = fopen(truth_set, "w");
	if (!fp) {
		printf("Could not create %s\n", truth_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  find ground truth results (using linear scan method)
	// -------------------------------------------------------------------------
	Result **result = new Result*[qn];
	for (int i = 0; i < qn; ++i) {
		result[i] = new Result[MAXK];
	}
	k_nn_search(n, qn, d, MAXK, p, data, query, result);

	fprintf(fp, "%d %d\n", qn, MAXK);
	for (int i = 0; i < qn; ++i) {
		for (int j = 0; j < MAXK; ++j) {
			fprintf(fp, "%d %f ", result[i][j].id_, result[i][j].key_);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < qn; ++i) {
		delete[] result[i]; result[i] = NULL;
	}
	delete[] result; result = NULL;

	gettimeofday(&g_end_time, NULL);
	float truth_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Ground Truth: %f Seconds\n\n", truth_time);

	return 0;
}

// -----------------------------------------------------------------------------
void k_nn_search(					// k-NN search
	int   n, 							// cardinality
	int   qn,							// query number
	int   d, 							// dimensionality
	int   k,							// top-k value
	float p,							// the p value of Lp norm, p in (0,2]
	const float **data,					// data objects
	const float **query,				// query objects
	Result **result)					// k-MIP results (return)
{
	using namespace std::chrono;
	// -------------------------------------------------------------------------
	//  k-NN search by linear scan
	// -------------------------------------------------------------------------
	MinK_List *list = new MinK_List(k);
	g_runtime2 = 0;
	for (int i = 0; i < qn; ++i) {
		float kdist = MAXREAL;
		list->reset();
		for (int j = 0; j < n; ++j) {
			float dist = calc_lp_dist(d, p, kdist, data[j], query[i]);
			kdist = list->insert(dist, j + 1);
		}

		for (int j = 0; j < k; ++j) {
			result[i][j].id_  = list->ith_id(j);
			result[i][j].key_ = list->ith_key(j);
		}
		
	}

	delete list; list = NULL;
}

// -----------------------------------------------------------------------------
void k_nn_search(					// k-NN search
	int   n, 							// cardinality
	int   qn,							// query number
	int   d, 							// dimensionality
	int   k,							// top-k value
	float p,							// the p value of Lp norm, p in (0,2]
	const char *data_file,					// data and query file on NVM
	Result **result)					// k-MIP results (return)
{
	using namespace std::chrono;
	// -------------------------------------------------------------------------
	//  k-NN search by linear scan
	// -------------------------------------------------------------------------
	MinK_List *list = new MinK_List(k);
	g_runtime2 = 0;

	// open file on NVM
	int is_pmem = 0;
	size_t mapped_len = 0;
	char* nvm_pool_ = (char*) pmem_map_file(data_file, DATA_POOL_SIZE, PMEM_FILE_CREATE, 0666, &mapped_len, &is_pmem);
	if (nvm_pool_ == NULL) {
		perror ("pmem_map_file");
		exit(1);
    }

    // printf("NVM mapping address: %p, size: %ld, is_pmem: %d\n", nvm_pool_, mapped_len, is_pmem);
	
	// query and data address
	char* query_addr = nvm_pool_;
	char* data_addr = nvm_pool_ + sizeof(float) * d * qn;
	float* query = new float[d];
	float* data = new float[d];
	
	for (int i = 0; i < qn; ++i) {
		// read query i
		g_start_time2 = high_resolution_clock::now();
		memcpy(query, query_addr + i * sizeof(float) * d, sizeof(float) * d);
		float kdist = MAXREAL;
		list->reset();
		for (int j = 0; j < n; ++j) {
			memcpy(data, data_addr + j * sizeof(float) * d, sizeof(float) * d);
			float dist = calc_lp_dist(d, p, kdist, data, query);
			kdist = list->insert(dist, j + 1);
		}

		for (int j = 0; j < k; ++j) {
			result[i][j].id_  = list->ith_id(j);
			result[i][j].key_ = list->ith_key(j);
		}
		g_end_time2 = high_resolution_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(g_end_time2 - g_start_time2);
		g_runtime2 += time_span.count();
	}
	g_runtime2 *= 1000;
	delete list; list = NULL;
	pmem_unmap(nvm_pool_, DATA_POOL_SIZE);
}

void write_nvm_data(
	int   n,
	int   qn,
	int   d,
	const float **data,					// data objects
	const float **query,				// query objects
	const char* file_path)
{
	// open file on NVM
	int is_pmem = 0;
	size_t mapped_len = 0;
	char* nvm_pool_ = (char*) pmem_map_file(file_path, DATA_POOL_SIZE, PMEM_FILE_CREATE, 0666, &mapped_len, &is_pmem);
	if (nvm_pool_ == NULL) {
		perror ("pmem_map_file");
		exit(1);
    }

    printf("NVM mapping address: %p, size: %ld, is_pmem: %d\n", nvm_pool_, mapped_len, is_pmem);

	// store query first
	char* tmp = nvm_pool_;
	for (int i = 0; i < qn; ++i) {
		memcpy(tmp + i * sizeof(float) * d, query[i], sizeof(float) * d);
	}
	// then store data
	tmp = nvm_pool_ + sizeof(float) * d * qn;
	for (int j = 0; j < n; ++j) {
		memcpy(tmp + j * sizeof(float) * d, data[j], sizeof(float) * d);
	}
	pmem_unmap(nvm_pool_, DATA_POOL_SIZE);
}

// -----------------------------------------------------------------------------
void test_rebuild_tree(					// k-NN search by qalsh
	int   n,							// number of data  objects
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	float zeta,							// symmetric factor of p-stable distr.
	float ratio,						// approximation ratio
	const float **data,					// data set
	const char*   data_file,            // data_file
	const char*   nvm_file,             // nvm file
	const char *out_path) 				// output path
{
	char output_set[200];
	sprintf(output_set, "%sqalsh_tree.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return;
	}
	QALSH_TREE *lsh = new QALSH_TREE(n, d, p, zeta, ratio, nvm_file, data_file, data);
	delete lsh; lsh = NULL;

	gettimeofday(&g_start_time, NULL);
	lsh = new QALSH_TREE(data_file, nvm_file);
	gettimeofday(&g_end_time, NULL);

	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Recover Time = %f Seconds\n", indexing_time);
	fprintf(fp, "Recover Time = %f Seconds\n", indexing_time);
}

// -----------------------------------------------------------------------------
int KVComp(						// compare function for qsort (ascending)
	const void *e1,						// 1st element
	const void *e2)						// 2nd element
{
	int ret = 0;
	KV *item1 = (KV*) e1;
	KV *item2 = (KV*) e2;

	if (item1->key_ < item2->key_) {
		ret = -1;
	} 
	else if (item1->key_ > item2->key_) {
		ret = 1;
	} 
	else {
		if (item1->id_ < item2->id_) ret = -1;
		else if (item1->id_ > item2->id_) ret = 1;
	}
	return ret;
}