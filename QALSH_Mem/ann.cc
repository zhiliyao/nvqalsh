#include "ann.h"
#include <vector>
using namespace std::chrono;

// -----------------------------------------------------------------------------
int linear_scan(					// k-NN search by linear scan
	int   n,							// number of data points
	int   qn,							// number of query points
	int   d,							// dimensionality
	float p,							// the p value of Lp norm, p in (0,2]
	const float **data,					// data set
	const float **query,				// query set
	const Result **R,					// truth set
	const char *out_path)				// output path
{
	char output_set[200];
	sprintf(output_set, "%slinear.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  k-NN search by linear scan
	// -------------------------------------------------------------------------
	Result **result = new Result*[qn];
	for (int i = 0; i < qn; ++i) result[i] = new Result[MAXK];

	printf("Top-k NN by Linear Scan:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		g_runtime2 = 0;
		g_start_time2 = high_resolution_clock::now();
		int top_k = TOPK[num];
		k_nn_search(n, qn, d, top_k, p, data, query, result);
		g_end_time2 = high_resolution_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(g_end_time2 - g_start_time2);
		g_runtime2 += time_span.count();
		g_runtime2 *= 1000;
		
		g_ratio  = 0.0f;
		g_recall = 0.0f;
		for (int i = 0; i < qn; ++i) {
			g_recall += calc_recall(top_k, R[i], (const Result *) result[i]);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += result[i][j].key_ / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}


		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime2 /= qn;

		printf("  %3d\t\t%.4f\t\t%.3lf\t\t%.2f\n", top_k, g_ratio, g_runtime2, 
			g_recall);
		fprintf(fp, "%d\t%f\t%lf\t%f\n", top_k, g_ratio, g_runtime2, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < qn; ++i) {
		delete[] result[i]; result[i] = NULL;
	}
	delete[] result; result = NULL;
	
	return 0;
}

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
	const char *out_path)				// output path
{	
	char output_set[200];
	sprintf(output_set, "%sqalsh.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	QALSH *lsh = new QALSH(n, d, p, zeta, ratio, data);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time = %f Seconds\n", indexing_time);
	printf("Memory = %f MB\n\n", g_memory / 1048576.0f);
	
	fprintf(fp, "index_time = %f Seconds\n", indexing_time);
	fprintf(fp, "memory     = %f MB\n\n", g_memory / 1048576.0f);

	// -------------------------------------------------------------------------
	//  c-k-ANN search
	// -------------------------------------------------------------------------
	printf("Top-k NN Search by QALSH:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall (%%)\n");
	for (int num = 0; num < MAX_ROUND; ++num) {
		int top_k = k;
		MinK_List *list = new MinK_List(top_k);

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		g_runtime2 = 0;
		g_start_time2 = high_resolution_clock::now();
		for (int i = 0; i < qn; ++i) {
			list->reset();
			
			lsh->knn(top_k, query[i], list);

			g_recall += calc_recall(top_k, (const Result*) R[i], list);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += list->ith_key(j) / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		g_end_time2 = high_resolution_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(g_end_time2 - g_start_time2);
		g_runtime2 = time_span.count(); // time in second
		delete list; list = NULL;

		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime2 = (g_runtime2 * 1000.0f) / qn;

		printf("  %3d\t\t%.4f\t\t%.3lf\t\t%.2f\n", top_k, g_ratio, g_runtime2, 
			g_recall);
		fprintf(fp, "%d\t%f\t%lf\t%f\n", top_k, g_ratio, g_runtime2, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

	return 0;
}

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
	const char *out_path) 				// output path
{	
	char output_set[200];
	sprintf(output_set, "%sqalsh_tree.out", out_path);

	FILE *fp = fopen(output_set, "a+");
	if (!fp) {
		printf("Could not create %s\n", output_set);
		return 1;
	}

	// -------------------------------------------------------------------------
	//  indexing
	// -------------------------------------------------------------------------
	gettimeofday(&g_start_time, NULL);
	QALSH_TREE *lsh = new QALSH_TREE(n, d, p, zeta, ratio, nvm_file, data_file, data);
	lsh->display();

	gettimeofday(&g_end_time, NULL);
	float indexing_time = g_end_time.tv_sec - g_start_time.tv_sec + 
		(g_end_time.tv_usec - g_start_time.tv_usec) / 1000000.0f;
	printf("Indexing Time = %f Seconds\n", indexing_time);
	printf("Memory = %f MB\n\n", g_memory / 1048576.0f);
	
	fprintf(fp, "index_time = %f Seconds\n", indexing_time);
	fprintf(fp, "memory     = %f MB\n\n", g_memory / 1048576.0f);

	// -------------------------------------------------------------------------
	//  c-k-ANN search
	// -------------------------------------------------------------------------
	printf("Top-k NN Search by QALSH_TREE:\n");
	printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall (%%)\n");
	std::vector<KV> cand; cand.resize(100);
	for (int num = 0; num < MAX_ROUND; ++num) {
		int top_k = k;

		g_ratio  = 0.0f;
		g_recall = 0.0f;
		g_runtime2 = 0;
		g_start_time2 = high_resolution_clock::now();
		for (int i = 0; i < qn; ++i) {
			lsh->knn(top_k, query[i], cand);

			g_recall += calc_recall(top_k, (const Result*) R[i], cand);

			float ratio = 0.0f;
			for (int j = 0; j < top_k; ++j) {
				ratio += cand[j].key_ / R[i][j].key_;
			}
			g_ratio += ratio / top_k;
		}
		g_end_time2 = high_resolution_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(g_end_time2 - g_start_time2);
		g_runtime2 += time_span.count(); // time in second
		g_ratio   = g_ratio / qn;
		g_recall  = g_recall / qn;
		g_runtime2 = (g_runtime2 * 1000) / qn; // time in ms

		printf("  %3d\t\t%.4f\t\t%.3lf\t\t%.2f\n", top_k, g_ratio, g_runtime2, 
			g_recall);
		fprintf(fp, "%d\t%f\t%lf\t%f\n", top_k, g_ratio, g_runtime2, g_recall);
	}
	printf("\n");
	fprintf(fp, "\n");
	fclose(fp);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete lsh; lsh = NULL;

	return 0;
}