#include "qalsh_tree.h"
#include <queue>

// some useful 
#define HASH_FUNC_SIZE ((int)ceil((double)(this->dim_) * (double)sizeof(float) / 256.0) * 256)

QALSH_TREE::QALSH_TREE(					// constructor
		int   n,						// cardinality
		int   d,						// dimensionality
		float p,						// l_p distance
		float zeta,						// symmetric factor of p-stable distr.
		float ratio,					// approximation ratio
		const char* nvm_file,           // nvm file
        const char* data_file,          // data file on disk
		const float** data)		    	// data objects
: n_pts_(n),
  dim_(d),
  p_(p),
  zeta_(zeta),
  ratio_(ratio),
  fp_(NULL),
  nvm_pool_(NULL)
{
	nvm_allocator = new NVMAllocator(nvm_file, POOL_SIZE);
	this->nvm_pool_ = (char *)nvmpool_alloc(256);
	worker_id= 0;
	the_thread_mempools.init(1, MEM_POOL_SIZE, 4096);

	init();
	
	// bulkload to build index array in memory
	bulkload(data);

	// open data file for future reding
	this->fp_ = fopen(data_file, "rb");
	if (this->fp_ == NULL) {
		perror ("open data file");
        exit(1);
	}
}

// -----------------------------------------------------------------------------
void QALSH_TREE::init()					// basic initialization
{
	// -------------------------------------------------------------------------
	//  init <w_> <m_> and <l_> (auto tuning-w)
	//  
	//  w0 ----- best w for L_{0.5} norm to minimize m (auto tuning-w)
	//  w1 ----- best w for L_{1.0} norm to minimize m (auto tuning-w)
	//  w2 ----- best w for L_{2.0} norm to minimize m (auto tuning-w)
	//  other w: use linear combination for interpolation
	// -------------------------------------------------------------------------
	float delta = 1.0f / E;
	float beta  = (float) CANDIDATES / (float) n_pts_;

	float w0 = (ratio_ - 1.0f) / log(sqrt(ratio_));
	float w1 = 2.0f * sqrt(ratio_);
	float w2 = sqrt((8.0f * SQR(ratio_) * log(ratio_)) / (SQR(ratio_) - 1.0f));
	float p1 = -1.0f, p2 = -1.0f;

	if (fabs(p_ - 0.5f) < FLOATZERO) {
		w_ = w0;
		p1 = calc_l0_prob(w_ / 2.0f);
		p2 = calc_l0_prob(w_ / (2.0f * ratio_));
	}
	else if (fabs(p_ - 1.0f) < FLOATZERO) {
		w_ = w1;
		p1 = calc_l1_prob(w_ / 2.0f);
		p2 = calc_l1_prob(w_ / (2.0f * ratio_));
	}
	else if (fabs(p_ - 2.0f) < FLOATZERO) {
		w_ = w2;
		p1 = calc_l2_prob(w_ / 2.0f);
		p2 = calc_l2_prob(w_ / (2.0f * ratio_));
	}
	else {
		if (fabs(p_-0.8f) < FLOATZERO) w_ = 2.503f;
		else if (fabs(p_-1.2f) < FLOATZERO) w_ = 3.151f;
		else if (fabs(p_-1.5f) < FLOATZERO) w_ = 3.465f;
		else w_ = (w2 - w1) * p_ + (2.0f * w1 - w2);

		new_stable_prob(p_, zeta_, ratio_, 1.0f, w_, 1000000, p1, p2);
	}

	float para1 = sqrt(log(2.0f / beta));
	float para2 = sqrt(log(1.0f / delta));
	float para3 = 2.0f * (p1 - p2) * (p1 - p2);
	float eta   = para1 / para2;
	float alpha = (eta * p1 + p2) / (1.0f + eta);

	m_ = (int) ceil((para1 + para2) * (para1 + para2) / para3);
	l_ = (int) ceil(alpha * m_);

	// -------------------------------------------------------------------------
	//  store meta-data in NVM
	// -------------------------------------------------------------------------
	// use first 256B to store parameter
	Root* root = (Root*) this->nvm_pool_;
	root->n_pts_ = this->n_pts_;
	root->dim_ = this->dim_;
	root->p_ = this->p_;
	root->zeta_ = this->zeta_;
	root->ratio_ = this->ratio_;
	root->w_ = this->w_;
	root->m_ = this->m_;
	root->l_ = this->l_;

	// -------------------------------------------------------------------------
	//  generate hash functions <a_>, and store them in NVM
	// -------------------------------------------------------------------------
	g_memory += SIZEFLOAT * m_ * dim_;
	a_ = new float*[m_];
	char* tmp_ptr;
	for (int i = 0; i < m_; ++i) {
		tmp_ptr = (char*)nvmpool_alloc(HASH_FUNC_SIZE);
		a_[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			if (fabs(p_-0.5f) < FLOATZERO) a_[i][j] = levy(1.0f, 0.0f);
			else if (fabs(p_-1.0f) < FLOATZERO) a_[i][j]= cauchy(1.0f, 0.0f);
			else if (fabs(p_-2.0f) < FLOATZERO) a_[i][j] = gaussian(0.0f, 1.0f);
			else a_[i][j] = p_stable(p_, zeta_, 1.0f, 0.0f);
		}
		memcpy(tmp_ptr, a_[i], sizeof(float) * dim_);
	}

	// generate {m_} hash table (B+-Trees)
	g_memory += sizeof(tree*) * this->m_;
	this->tables_ = new tree*[this->m_];
	char* begin = (char *)nvmpool_alloc(256);
	char* tmp = begin;
	for (int i = 0; i < this->m_; i++) {
		if ((tmp - begin) >= 256) {
			begin = (char *)nvmpool_alloc(256);
			tmp = begin;
		}
		this->tables_[i] = initTree(tmp, false);
		uint64_t* pp = (uint64_t*) tmp;
		uint64_t* pp2 = (uint64_t*) (tmp + sizeof(uint64_t));
		tmp += 16;
	}
}

QALSH_TREE::QALSH_TREE(
		const char* data_file,
		const char* nvm_file)
: fp_(NULL),
  nvm_pool_(NULL)
{
	nvm_allocator = new NVMAllocator(nvm_file, POOL_SIZE);
	this->nvm_pool_ = (char *)nvmpool_alloc(256);
	worker_id= 0;
	the_thread_mempools.init(1, MEM_POOL_SIZE, 4096);

	// restore the parameters
	Root* root = (Root*) this->nvm_pool_;
	this->n_pts_ = root->n_pts_;
	this->dim_ = root->dim_;
	this->p_ = root->p_;
	this->zeta_ = root->zeta_;
	this->ratio_ = root->ratio_;
	this->w_ = root->w_;
	this->m_ = root->m_;
	this->l_ = root->l_;

	// restore hash functions
	g_memory += SIZEFLOAT * m_ * dim_;
	a_ = new float*[m_];
	char* tmp_ptr;
	for (int i = 0; i < m_; ++i) {
		a_[i] = new float[dim_];
		tmp_ptr = (char*)nvmpool_alloc(HASH_FUNC_SIZE);
		memcpy(a_[i], tmp_ptr, sizeof(float) * dim_);
	}

	// generate {m_} hash table (B+-Trees)
	g_memory += sizeof(tree*) * this->m_;
	this->tables_ = new tree*[this->m_];
	char* begin = (char *)nvmpool_alloc(256);
	char* tmp = begin;
	for (int i = 0; i < this->m_; i++) {
		if ((tmp - begin) >= 256) {
			begin = (char *)nvmpool_alloc(256);
			tmp = begin;
		}
		uint64_t* pp = (uint64_t*) tmp;
		uint64_t* pp2 = (uint64_t*) (tmp + sizeof(uint64_t));
		this->tables_[i] = initTree(tmp, true);
		this->tables_[i]->reBulkload(g_memory);
		tmp += 16;
	}

	// open data file for future reding
	this->fp_ = fopen(data_file, "rb");
	if (this->fp_ == NULL) {
		perror ("open data file");
        exit(1);
	}
}

void QALSH_TREE::bulkload(const float** data) {

	// -------------------------------------------------------------------------
	//  bulkloading
	// -------------------------------------------------------------------------
	KV* table = new KV[this->n_pts_];
	for (int i = 0; i < m_; ++i) {
		for (int j = 0; j < this->n_pts_; ++j) {
			table[j].id_  = j;
			table[j].key_ = calc_hash_value(i, data[j]);
		}
		qsort(table, this->n_pts_, sizeof(KV), KVComp);
		this->tables_[i]->bulkload(this->n_pts_, table, 1.0, g_memory);
	}
	delete[] table;
}

QALSH_TREE::~QALSH_TREE()						// destructor
{
	for (int i = 0; i < m_; ++i) {
		delete[] a_[i]; a_[i] = NULL;
		delete tables_[i]; tables_[i] = NULL;
	}
	delete[] a_; a_ = NULL;
	delete[] tables_; tables_ = NULL;

	g_memory -= SIZEFLOAT * m_ * dim_;
	g_memory -= this->m_ * sizeof(lbtree*);

	// unmap NVM
	delete nvm_allocator;

	// close data file
	fclose(this->fp_);
}

// -----------------------------------------------------------------------------
inline float QALSH_TREE::calc_l0_prob(	// calc prob of L1/2 dist
	float x)							// x = w / (2.0 * r)
{
	return new_levy_prob(x);
}

// -----------------------------------------------------------------------------
inline float QALSH_TREE::calc_l1_prob(	// calc prob of L1 dist
	float x)							// x = w / (2.0 * r)
{
	return new_cauchy_prob(x);
}

// -----------------------------------------------------------------------------
inline float QALSH_TREE::calc_l2_prob(	// calc prob of L2 dist
	float x)							// x = w / (2.0 * r)
{
	return new_gaussian_prob(x);
}

// -----------------------------------------------------------------------------
float QALSH_TREE::calc_hash_value( 		// calc hash value
	int   tid,							// hash table id
	const float *data)					// one data/query object
{
	return calc_inner_product(dim_, a_[tid], data);
}

// -----------------------------------------------------------------------------
void QALSH_TREE::display()				// display parameters
{
	printf("Parameters of QALSH_TREE:\n");
	printf("    n     = %d\n",   n_pts_);
	printf("    d     = %d\n",   dim_);
	printf("    p     = %.1f\n", p_);
	printf("    zeta  = %.1f\n", zeta_);
	printf("    ratio = %.1f\n", ratio_);
	printf("    w     = %f\n",   w_);
	printf("    m     = %d\n",   m_);
	printf("    l     = %d\n",   l_);
	printf("\n");
}

int QALSH_TREE::knn(						// k-NN search
    int   top_k,					// top-k value
    const float *query,				// input query object
    std::vector<KV> &cand)		// k-NN results (return)
{
    int   *freq    = new int[n_pts_];
    // left {m_} bleaf
	bleaf* *lpos   = new bleaf*[m_];
    // right {m_} belaf
	bleaf* *rpos   = new bleaf*[m_];
    // leaf/right projection distance
    // float* ldist   = new float[m_];
    // float* rdist   = new float[m_];
	bool  *checked = new bool[n_pts_];
	bool  *flag    = new bool[m_];
	float *q_val   = new float[m_];
	float *data    = new float[dim_];

	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	memset(freq,    0,     n_pts_ * SIZEINT);
	memset(checked, false, n_pts_ * SIZEBOOL);

    int tmp;
	for (int i = 0; i < m_; ++i) {
		q_val[i] = calc_hash_value(i, query);
        lpos[i] = (bleaf*)tables_[i]->lookup(q_val[i], &tmp);
        rpos[i] = lpos[i]->getNextLeaf();
	}

	// -------------------------------------------------------------------------
	//  c-k-ANN search
	// -------------------------------------------------------------------------
	int candidates = CANDIDATES + top_k - 1; // candidate size
	int cand_cnt   = 0;				// candidate counter
	
	float kdist  = MAXREAL;
	float radius = 1.0f;			// search radius
	float width  = radius * w_ / 2.0f; // bucket width

	std::priority_queue<KV> myqueue;
	int total_io = 0;
	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialize the stop condition for current round
		// ---------------------------------------------------------------------
		int num_flag = 0;
		memset(flag, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: (R,c)-NN search
		// ---------------------------------------------------------------------
		while (num_flag < m_) {
			for (int j = 0; j < m_; ++j) {
				if (!flag[j]) continue;

				// table = tables_[j];
				float q_v = q_val[j], ldist = -1.0f, rdist = -1.0f;
				// -------------------------------------------------------------
				//  step 2.1: scan the left part of hash table
				// -------------------------------------------------------------
				int cnt = 0;
				bleaf* pos = lpos[j];
				while (cnt < SCAN_SIZE_TREE) {
					ldist = MAXREAL;
					if (pos != NULL) {
						ldist = fabs(q_v - pos->getMinv());
					}
					else break;
					if (ldist > width) break;

                    uint32_t bitmap = pos->bitmap;
                    while (bitmap) {
                        int in_pos = __builtin_ffs(bitmap)-1;
                        int id = pos->ch(in_pos);
                        if (++freq[id] >= l_ && !checked[id]) {
                            checked[id] = true;
                            getDataById(id, data);
                            total_io++;
                            float dist = calc_lp_dist(dim_, p_, kdist, data, query);
                            if (myqueue.size() < top_k) {
                                myqueue.push({dist, id});
                            } else {
                                if (dist < myqueue.top().key_) {
                                    myqueue.pop();
                                    myqueue.push({dist, id});
                                    kdist = myqueue.top().key_;
                                }
                            }

                            if (++cand_cnt >= candidates) break;
                        }
						bitmap &= ~(0x1<<in_pos);
                    }
                    if (cand_cnt >= candidates) break;

					pos = pos->getLastLeaf(); ++cnt;
				}
				if (cand_cnt >= candidates) break;
				lpos[j] = pos;

				// -------------------------------------------------------------
				//  step 2.2: scan the right part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				pos = rpos[j];
				while (cnt < SCAN_SIZE_TREE) {
					rdist = MAXREAL;
					if (pos != NULL) {
						rdist = fabs(q_v - pos->getMaxv());
					}
					else break;
					if (rdist > width) break;

					uint32_t bitmap = pos->bitmap;
                    while (bitmap) {
                        int in_pos = __builtin_ffs(bitmap)-1;
                        int id = pos->ch(in_pos);
                        if (++freq[id] >= l_ && !checked[id]) {
                            checked[id] = true;
                            getDataById(id, data);
                            total_io++;
                            float dist = calc_lp_dist(dim_, p_, kdist, data, query);
                            if (myqueue.size() < top_k) {
                                myqueue.push({dist, id});
                            } else {
                                if (dist < myqueue.top().key_) {
                                    myqueue.pop();
                                    myqueue.push({dist, id});
                                    kdist = myqueue.top().key_;
                                }
                            }

                            if (++cand_cnt >= candidates) break;
                        }
						bitmap &= ~(0x1<<in_pos);
                    }
                    if (cand_cnt >= candidates) break;

					pos = pos->getNextLeaf(); ++cnt;
				}
				if (cand_cnt >= candidates) break;
				rpos[j] = pos;

				// -------------------------------------------------------------
				//  step 2.3: check whether this width is finished scanned
				// -------------------------------------------------------------
				if (ldist > width && rdist > width) {
					flag[j] = false;
					if (++num_flag >= m_) break;
				}
			}
			if (num_flag >= m_ || cand_cnt >= candidates) break;
		}
		// ---------------------------------------------------------------------
		//  step 3: stop conditions t1 and t2
		// ---------------------------------------------------------------------
		if (kdist < ratio_ * radius && cand_cnt >= top_k) break;
		if (cand_cnt >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: auto-update radius
		// ---------------------------------------------------------------------
		radius = radius * ratio_;
		width  = radius * w_ / 2.0f;
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] freq;    freq    = NULL;
	delete[] lpos;    lpos    = NULL;
	delete[] rpos;    rpos    = NULL;
	delete[] checked; checked = NULL;
	delete[] flag;    flag    = NULL;
	delete[] q_val;   q_val   = NULL;
	delete[] data;    data    = NULL;

	int ii = top_k - 1;
	while (!myqueue.empty()) {
		cand[ii--] = myqueue.top();
		myqueue.pop();
	}

	return total_io;
	return 0;
}

void QALSH_TREE::getDataById(int id, float* data) {
	if (fseek(this->fp_, id * sizeof(float) * this->dim_, SEEK_SET) != 0) {
		perror("fseek");
		exit(-1);
	}
	if (fread(data, sizeof(float), this->dim_, this->fp_) != this->dim_) {
		perror("read data");
		exit(-1);
	}
}

void QALSH_TREE::check() {
	printf("n: %d, d: %d, m: %d\n", n_pts_, dim_, m_);
	int tmp = -1;
	bleaf* bl1 = (bleaf*)tables_[0]->lookup(50, &tmp);
	printf("%d %d\n", ((char*)bl1 - nvm_pool_), tmp);
	bleaf* bl2 = (bleaf*)tables_[m_-1]->lookup(50, &tmp);
	printf("%d %d\n", ((char*)bl2 - nvm_pool_), tmp);
	printf("hf: %f, %f\n", a_[0][0], a_[m_-1][dim_-1]);
	printf("%f\n", (g_memory / 1024.0 / 1024.0));
}