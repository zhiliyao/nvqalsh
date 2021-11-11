/**
 * @file indirect-arronly.cc
 * @author  Shimin Chen <shimin.chen@gmail.com>, Jihang Liu, Leying Chen
 * @version 1.0
 *
 * @section LICENSE
 *
 * TBD
 *
 * @section DESCRIPTION
 *
 *
 * The class implements a Btree with indirect arrays only.  Each node contains
 * an indirect array.  The design aims to have good search performance and
 * efficient solution for update consistency on NVM memory.  However, the node
 * size is limited to up to 128B.
 */

#include "lbtree.h"

/* ----------------------------------------------------------------- *
 useful structure
 * ----------------------------------------------------------------- */
static int last_slot_in_line[LEAF_KEY_NUM];

static void initUseful(void)
{
    // line 0
    last_slot_in_line[0]= 2;
    last_slot_in_line[1]= 2;
    last_slot_in_line[2]= 2;

    // line 1
    last_slot_in_line[3]= 6;
    last_slot_in_line[4]= 6;
    last_slot_in_line[5]= 6;
    last_slot_in_line[6]= 6;

    // line 2
    last_slot_in_line[7]= 10;
    last_slot_in_line[8]= 10;
    last_slot_in_line[9]= 10;
    last_slot_in_line[10]=10;

    // line 3
    last_slot_in_line[11]=13;
    last_slot_in_line[12]=13;
    last_slot_in_line[13]=13;
}

/* ----------------------------------------------------------------- *
 bulk load
 * ----------------------------------------------------------------- */

/* generate a btree on the input keys.
 After generating the btree, level is returned.
 leaf and non-leaf nodes will be bfill full.
 bfill should be between 0 and 1.
 */


/**
 * build a subtree containing: input[start_key, start_key+num_key-1]
 * stop at target_level  (leaf is at level 0)
 *
 * @param input          keyInput instance 
 * @param start_key      keyInput index for the first key
 * @param num_key        number of keys in the subtree
 * @param bfill          filling factor in (0.0,1.0] 
 * @param target_level   stop buidling the subtree at target level
 * @param pfirst         return pointers to the first node at each level
 * @param n_nodes        return number of nodes at each level 
 *
 * @retval the top level of the subtree (can be smaller than target_level
 * if the root level is smaller target_level)
 */
int lbtree::bulkloadSubtree(
    KV *input, int start_key, int num_key, 
    float bfill, int target_level,
    Pointer8B pfirst[], int n_nodes[], uint64_t& g_memory)
{
    // We assume that a tree cannot be higher than 32 levels
    int ncur[32];     // current node at every level
    int top_level;    // top_level is the top level that this method builds

    assert(start_key>=0 && num_key > 0 && bfill>0.0 && bfill<=1.0
           && target_level>=0);


    // 1. compute leaf and nonleaf number of keys
    int leaf_fill_num= (int)((float)LEAF_KEY_NUM * bfill);
    leaf_fill_num= Max(leaf_fill_num, 1);

    int nonleaf_fill_num= (int)((float)NON_LEAF_KEY_NUM * bfill);
    nonleaf_fill_num= Max(nonleaf_fill_num, 1);


    // 2. compute number of nodes
    n_nodes[0]= ceiling(num_key, leaf_fill_num);
    top_level= 0;

    for (int i=1; n_nodes[i-1]>1 && i<=target_level; i++) {
       n_nodes[i]= ceiling(n_nodes[i-1], nonleaf_fill_num+1);
       top_level= i;
    } // end of for each nonleaf level


    // 3. allocate nodes
    pfirst[0]= nvmpool_alloc(sizeof(bleaf) * n_nodes[0]);
    for (int i=1; i<=top_level; i++) {
        g_memory += sizeof(bnode) * n_nodes[i];
       pfirst[i]= mempool_alloc(sizeof(bnode) * n_nodes[i]);
    }
       // nvmpool_alloc/mempool_alloc call exit if out of memory


    // 4. populate nodes
    for (int ll=1; ll<=top_level; ll++) {
        ncur[ll]= 0;
        bnode *np= (bnode *)(pfirst[ll]);
        np->lock()= 0; np->num()= -1;
    }

    bleaf * leaf= pfirst[0];
    int nodenum= n_nodes[0];

    bleafMeta leaf_meta;
    leaf_meta.v.bitmap= ( ((1<<leaf_fill_num)-1)
                         <<(LEAF_KEY_NUM-leaf_fill_num));
    leaf_meta.v.lock= 0;
    leaf_meta.v.alt= 0;

    int key_id= start_key;

    bleaf* last_leaf = NULL;
    for (int i=0; i<nodenum; i++) {
        bleaf *lp= &(leaf[i]);

        // compute number of keys in this leaf node
        int fillnum= leaf_fill_num; // in most cases
        if (i == nodenum-1) {
           fillnum= num_key - (nodenum-1)*leaf_fill_num;
           assert(fillnum>=1 && fillnum<=leaf_fill_num);

           leaf_meta.v.bitmap= ( ((1<<fillnum)-1)
                                <<(LEAF_KEY_NUM-fillnum));
        }
        leaf_meta.v.min_idx = LEAF_KEY_NUM - fillnum;
        // lbtree tends to leave the first line empty
        for (int j=LEAF_KEY_NUM-fillnum; j<LEAF_KEY_NUM; j++) {
            key_id ++;

            // entry
            lp->k(j) = input[key_id].key_;
            lp->ch(j) = input[key_id].id_;
            

        } // for each key in this leaf node
        leaf_meta.v.max_idx = LEAF_KEY_NUM - 1;
        // sibling pointer
        lp->setNextLeaf(((i<nodenum-1) ? &(leaf[i+1]) : NULL));
        lp->setlastLeaf(last_leaf);

        // 2x8B meta
        lp->setBothHalfWord(&leaf_meta);


        // populate nonleaf node
        Pointer8B child= lp;
        key_type  left_key= lp->k(LEAF_KEY_NUM-fillnum);

        // append (left_key, child) to level ll node
        // child is the level ll-1 node to be appended.
        // left_key is the left-most key in the subtree of child.
        for (int ll=1; ll<=top_level; ll++) {
           bnode *np= ((bnode *)(pfirst[ll])) + ncur[ll];

           // if the node has >=1 child
           if (np->num() >= 0) {
               int kk= np->num()+1;
               np->ch(kk)= child; np->k(kk)= left_key; 
               np->num()= kk;

               if ((kk==nonleaf_fill_num)&&(ncur[ll]<n_nodes[ll]-1)) { 
                   ncur[ll] ++; np ++;
                   np->lock()= 0; np->num()= -1;
               }
               break;
           }

           // new node
           np->ch(0)= child; np->num() = 0;

           // append this new node to parent
           child= np;

        } // for each nonleaf level
        last_leaf = lp;
    } // end of foreach leaf node

    // 5. return
    return top_level;
}

/**
 * Obtain the node pointers and left keys for a given level
 *
 * @param pnode          subtree root
 * @param pnode_level    subtree root level
 * @param left_key       left key of the subtree
 * @param target_level   the level to get nodes and keys
 * @param ptrs           child node pointers (output)
 * @param keys           left keys of subtrees rooted at child nodes (output)
 * @param num_nodes      number of child nodes (output)
 * @param free_above_level_nodes  if nodes above target_level should be freed
 *
 */
void lbtree::getKeyPtrLevel(
Pointer8B pnode, int pnode_level, key_type left_key,
int target_level, Pointer8B ptrs[], key_type keys[], int &num_nodes,
bool free_above_level_nodes)
{
    // already at this target_level
    if (pnode_level == target_level) {
        ptrs[num_nodes]= pnode;
        keys[num_nodes]= left_key;
        num_nodes ++;
        return;
    }

    // pnode_level > target_level
    else if (pnode_level > target_level) {
        bnode *p= pnode;
        getKeyPtrLevel(p->ch(0), pnode_level-1, left_key,
                       target_level, ptrs, keys, num_nodes,
                       free_above_level_nodes);
        for (int i=1; i<=p->num(); i++) {
            getKeyPtrLevel(p->ch(i), pnode_level-1, p->k(i),
                       target_level, ptrs, keys, num_nodes,
                       free_above_level_nodes);
        }

        if (free_above_level_nodes) mempool_free_node(p);
    }
}


typedef struct BldThArgs {

    key_type start_key; // input
    key_type num_key;   // input

    int top_level;  // output
    int n_nodes[32];  // output
    Pointer8B pfirst[32];  // output

} BldThArgs;

int lbtree::bulkload (int keynum, KV *input, float bfill, uint64_t& g_memory)
{
    // 1. allocate BldThArgs[]
    int num_threads = 1;

    BldThArgs *bta= new BldThArgs[num_threads];
    if (!bta) {perror("malloc"); exit(1);}


    // 2. one thread?
    if (num_threads == 1) {
        bta[0].top_level= bulkloadSubtree(
                               input, 0, keynum, bfill, 31,
                               bta[0].pfirst, bta[0].n_nodes, g_memory);
        tree_meta->root_level=  bta[0].top_level;
        tree_meta->tree_root=   bta[0].pfirst[tree_meta->root_level];
        tree_meta->setFirstLeaf(bta[0].pfirst[0]);
        tree_meta->setLeafNum(bta[0].n_nodes[0]);
        // if this assertion is false, then the tree has > 31 levels
        assert(bta[0].n_nodes[bta[0].top_level] == 1);

        delete[] bta;
        return tree_meta->root_level;
    }
}

int lbtree::reBulkload(uint64_t& g_memory) {
    uint64_t ppointer = *(tree_meta->first_leaf);
    uint64_t leaf_num = *(tree_meta->leaf_num);
    Pointer8B pfirst[32];
    int n_nodes[32];
    int ncur[32];
    int top_level;
    int target_level = 32;
    
    int nonleaf_fill_num= (int)((float)NON_LEAF_KEY_NUM * 1);
    nonleaf_fill_num= Max(nonleaf_fill_num, 1);

    if (ppointer == 0) {
        return -1;
    }
    bleaf** leafArray = (bleaf**)nvm_allocator->getVirtualAddr(ppointer);

    // compute number of nodes
    n_nodes[0] = leaf_num;
    top_level = 0;
    for (int i = 1; n_nodes[i-1] > 1 && i <= target_level; i++) {
        n_nodes[i] = ceiling(n_nodes[i-1], nonleaf_fill_num+1);
        top_level = i;
    }

    pfirst[0] = leafArray;
    // allocate inner nodes
    for (int i = 1; i <= top_level; i++) {
        g_memory += sizeof(bnode) * n_nodes[i];
        pfirst[i] = mempool_alloc(sizeof(bnode) * n_nodes[i]);
    }
    // populate nodes
    for (int ll = 1; ll <= top_level; ll++) {
        ncur[ll] = 0;
        bnode *np = (bnode*)(pfirst[ll]);
        np->lock() = 0; np->num() = -1;
    }
    
    bleaf* leaf = pfirst[0];
    // iterate the leaves nodes
    for (int i = 0; i < leaf_num; i++) {
        bleaf *lp = &(leaf[i]);

        // populate the nonleaf node
        Pointer8B child = lp;
        int first_pos= bitScan(lp->bitmap)-1;
        key_type left_key = lp->k(first_pos);

        // append (left_key, child) to level ll node
        // child is the level ll-1 node to be appended.
        // left_key is the left-most key in the subtree of child.
        for (int ll=1; ll<=top_level; ll++) {
           bnode *np= ((bnode *)(pfirst[ll])) + ncur[ll];

           // if the node has >=1 child
           if (np->num() >= 0) {
               int kk= np->num()+1;
               np->ch(kk)= child; np->k(kk)= left_key; 
               np->num()= kk;

               if ((kk==nonleaf_fill_num)&&(ncur[ll]<n_nodes[ll]-1)) { 
                   ncur[ll] ++; np ++;
                   np->lock()= 0; np->num()= -1;
               }
               break;
           }

           // new node
           np->ch(0)= child; np->num() = 0;

           // append this new node to parent
           child= np;

        } // for each nonleaf level
    }
    // set the metadata of tree
    tree_meta->root_level=  top_level;
    tree_meta->tree_root=   pfirst[top_level];
    return top_level;
}

/* ----------------------------------------------------------------- *
 look up
 * ----------------------------------------------------------------- */

/* leaf is level 0, root is level depth-1 */

void * lbtree::lookup (key_type key, int *pos)
{
    bnode *p;
    bleaf *lp;
    int i,t,m,b;
    key_type r;
    
    // unsigned char key_hash= hashcode1B(key);
    int ret_pos;
    
Again1:

    // 2. search nonleaf nodes
    p = tree_meta->tree_root;
    
    for (i=tree_meta->root_level; i>0; i--) {
        
        // prefetch the entire node
        NODE_PREF(p);
        
        // binary search to narrow down to at most 8 entries
        b=1; t=p->num();
        while (b+7<=t) {
            m=(b+t) >>1;
            r= key - p->k(m);
            if (r>0) b=m+1;
            else if (r<0) t = m-1;
            else {p=p->ch(m); goto inner_done;}
        }
        
        // sequential search (which is slightly faster now)
        for (; b<=t; b++)
            if (key < p->k(b)) break;
        p = p->ch(b-1);
        
    inner_done: ;
    }
    
    // 3. search leaf node
    lp= (bleaf *)p;

    return (void *)lp;
}

/* ------------------------------------- *
   quick sort the keys in leaf node
 * ------------------------------------- */


// pos[] will contain sorted positions
void lbtree::qsortBleaf(bleaf *p, int start, int end, int pos[])
{
    if (start >= end) return;

    int pos_start= pos[start];
    key_type key= p->k(pos_start);  // pivot
    int l, r;
        
    l= start;  r=end;
    while (l<r) {
        while ((l<r) && (p->k(pos[r])>key)) r--;
        if (l<r) {
            pos[l]= pos[r];
            l++;
        }
        while ((l<r) && (p->k(pos[l])<=key)) l++;
        if (l<r) {
            pos[r]= pos[l];
            r--;
        }
    }
    pos[l]= pos_start;
    qsortBleaf(p, start, l-1, pos);
    qsortBleaf(p, l+1, end, pos);
}

/* ----------------------------------------------------------------- *
 print
 * ----------------------------------------------------------------- */
void lbtree::print (Pointer8B pnode, int level)
{
    if (level > 0) {
        bnode *p= pnode;

        printf("%*cnonleaf lev=%d num=%d\n", 10+level*4, '+', level, p->num());

        print (p->ch(0), level-1);
        for (int i=1; i<=p->num(); i++) {
            printf ("%*c%lld\n", 10+level*4, '+', p->k(i));
            print (p->ch(i), level-1);
        }
    }
    else {
        bleaf * lp = pnode;
        printf("Min Key Index = %d\n", lp->min_idx);
        printf("Max Key Index = %d\n", lp->max_idx);
        uint32_t bmp= lp->bitmap;
        for (int i=0; i<LEAF_KEY_NUM; i++) {
            if (bmp & (1<<i)) {
                printf ("[%2d] key=%f\n", i, lp->k(i));
            }
        }

        bleaf * pnext= lp->nextSibling();
        if (pnext != NULL) {
            int first_pos= bitScan(pnext->bitmap)-1;
            printf ("->(%f)\n", pnext->k(first_pos));
        }
        else
            printf ("->(null)\n");
    }
}

/**
 * get min and max key in the given leaf p
 */
void lbtree::getMinMaxKey (bleaf *p, key_type &min_key, key_type &max_key)
{
    unsigned short bmp= p->bitmap;
    max_key= MIN_KEY;
    min_key= MAX_KEY;
    
    for (int i=0; i<LEAF_KEY_NUM; i++) {
        if (bmp & (1<<i)) {
            if (p->k(i) > max_key) max_key = p->k(i);
            if (p->k(i) < min_key) min_key = p->k(i);
        }
    }
}

void lbtree::checkFirstLeaf(void)
{
    // get left-most leaf node
    bnode *p= tree_meta->tree_root;
    for (int i= tree_meta->root_level; i>0; i--) p= p->ch(0);

    if ((bleaf *)p != nvm_allocator->getVirtualAddr(*tree_meta->first_leaf)) {
       printf("first leaf %p != %p\n", nvm_allocator->getVirtualAddr(*tree_meta->first_leaf), p);
       exit(1);
    }
}

/**
 * recursively check the subtree rooted at pnode
 *
 * If it encounters an error, the method will print an error message and exit.
 *
 * @param pnode   the subtree root
 * @param level   the level of pnode
 * @param start   return the start key of this subtree
 * @param end     return the end key of this subtree
 * @param ptr     ptr is the leaf before this subtree. 
 *                Upon return, ptr is the last leaf of this subtree.
 */
void lbtree::check (Pointer8B pnode, int level, key_type &start, key_type &end, bleaf * &ptr)
{
    if (pnode.isNull()) {
        printf ("level %d: null child pointer\n", level + 1);
        exit (1);
    }
    
    if (level == 0) { // leaf node
        bleaf *lp = pnode;

        if (((unsigned long long)lp)%256 != 0) {
            printf ("leaf(%p): not aligned at 256B\n", lp); exit (1);
        }
        
        // check number of keys
        if (lp->num() < 1) {// empty node!
            printf ("leaf(%p): empty\n", lp); exit (1);
        }

        // get min max
        getMinMaxKey(lp, start, end);

        // check lock bit
        if (lp->lock != 0) {
            printf ("leaf(%lld): lock bit == 1\n", start); 
            exit(1);
        }

        // check sibling pointer
        if ((ptr) && (ptr->nextSibling()!=lp)) {
            printf ("leaf(%lld): sibling broken from previous node\n", start); 
            fflush(stdout);

            /* output more info */
            bleaf *pp= (bleaf *)ptr;
            key_type ss, ee;
            getMinMaxKey(pp, ss, ee);
            printf("previous(%lld - %lld) -> ", ss, ee);

            pp= pp->nextSibling();
            if (pp == NULL) {
                printf("nil\n");
            }
            else {
                getMinMaxKey(pp, ss, ee);
                printf("(%lld - %lld)\n", ss, ee);
            }

            exit(1);
        }

        ptr= lp;
    }

    else { // nonleaf node
        key_type curstart, curend;
        int i;
        bleaf *curptr;
        
        bnode * p = pnode;

        if (((unsigned long long)p)%64!= 0) {
            printf ("nonleaf level %d(%p): not aligned at 64B\n", level, p); exit (1);
        }
        
        // check num of keys
        if (p->num()<0) {
            printf ("nonleaf level %d(%p): num<0\n", level, p); exit (1);
        }

        // check child 0
        curptr = ptr;
        check (p->ch(0), level-1, curstart, curend, curptr);
        start = curstart;
        if (p->num()>=1 && curend >= p->k(1)) {
            printf ("nonleaf level %d(%lld): key order wrong at child 0\n", level, p->k(1));
            exit (1);
        }
        
        // check child 1..num-1
        for (i=1; i<p->num(); i++) {
            check (p->ch(i), level-1, curstart, curend, curptr);
            if (!(p->k(i)<=curstart && curend<p->k(i+1)) ) {
                printf ("nonleaf level %d(%lld): key order wrong at child %d(%lld)\n",
                        level, p->k(1), i, p->k(i));
                exit (1);
            }
        }

        // check child num (when num>=1)
        if (i == p->num()) {
           check (p->ch(i), level-1, curstart, curend, curptr);
           if (curstart < p->k(i)) {
               printf ("nonleaf level %d(%lld): key order wrong at last child %d(%lld)\n", 
                        level, p->k(1), i, p->k(i));
               exit (1);
           }
        }
        end = curend;

        // check lock bit
        if (p->lock() != 0) {
            printf ("nonleaf level %d(%lld): lock bit is set\n", level, p->k(1)); 
            exit(1);
        }
        
        ptr = curptr;
    }
}


/* ------------------------------------------------------------------------- */
/*                              driver                                       */
/* ------------------------------------------------------------------------- */
tree * initTree(void *nvm_addr, bool recover)
{
    tree *mytree = new lbtree(nvm_addr, recover);
    return mytree;
}