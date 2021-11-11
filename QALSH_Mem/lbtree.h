/**
 * @file lbtree.h
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
 * The class implements the LB+-Tree. 
 *
 * Non-leaf nodes are in DRAM.  They are normal B+-Tree nodes.
 * Leaf nodes are in NVM.
 */

#ifndef _LBTREE_H
#define _LBTREE_H
/* ---------------------------------------------------------------------- */

#include "tree.h"

/* ---------------------------------------------------------------------- */

/* In a non-leaf, there are NON_LEAF_KEY_NUM keys and NON_LEAF_KEY_NUM+1
 * child pointers.
 */
#define NON_LEAF_KEY_NUM    (NONLEAF_SIZE/(KEY_SIZE+POINTER_SIZE)-1)

/* In a leaf, there are 8B header, 29x8B entries, 2x8B sibling pointers.
 */
#if LEAF_SIZE != 256
#error "LB+-Tree requires leaf node size to be 256B."
#endif

#define LEAF_KEY_NUM        (29) 

/* ---------------------------------------------------------------------- */
/**
 * Pointer8B defines a class that can be assigned to either bnode or bleaf.
 */
class Pointer8B {
public:
    unsigned long long  value;  /* 8B to contain a pointer */
    
public:
    Pointer8B(){}

    Pointer8B(const void *ptr)
    { value= (unsigned long long)ptr; }

    Pointer8B(const Pointer8B & p)
    { value= p.value; }

    Pointer8B & operator= (const void * ptr)
    {
        value= (unsigned long long)ptr;
        return *this;
    }
    Pointer8B & operator= (const Pointer8B & p)
    {
        value= p.value;
        return *this;
    }

    bool operator== (const void * ptr){
	bool result = (value==(unsigned long long)ptr);
        return result;
    }
    bool operator== (const Pointer8B & p){
	bool result = (value==p.value);
        return result;
    }
    
    
    operator void*() { return (void *)value; }
    operator char*() { return (char *)value; }
    operator struct bnode *() { return (struct bnode *)value;}
    operator struct bleaf *() { return (struct bleaf *)value;}
    operator unsigned long long () { return value;}
    
    bool isNull(void) {return (value==0);}
 
    void print(void) { printf("%llx\n", value); }
    
}; // Pointer8B

/**
 *  An IdxEntry consists of a key and a pointer.
 */
typedef struct IdxEntry{
    key_type   k;
    Pointer8B  ch;
} IdxEntry;

/**
 *  bnodeMeta: the 8B meta data in Non-leaf node
 */
typedef struct bnodeMeta {/* 8B */
    int  lock;    /* lock bit for concurrency control */
    int  num;     /* number of keys */
} bnodeMeta;

/**
 * bnode: non-leaf node
 *
 *   metadata (i.e. k(0))
 *
 *      k(1) .. k(NON_LEAF_KEY_NUM)
 *
 *   ch(0), ch(1) .. ch(NON_LEAF_KEY_NUM)
 */
class bnode{
public:
    IdxEntry   ent[NON_LEAF_KEY_NUM + 1];
public:
    key_type  & k(int idx)  { return ent[idx].k; }
    Pointer8B & ch(int idx) { return ent[idx].ch; }

    char * chEndAddr(int idx) {
       return (char *)&(ent[idx].ch)+sizeof(Pointer8B)-1;
    }
    
    int & num(void)  { return ((bnodeMeta *)&(ent[0].k))->num;}
    int & lock(void) { return ((bnodeMeta *)&(ent[0].k))->lock;}
}; // bnode

typedef union bleafMeta {
    uint32_t halfWord4B[2];
    struct {
        uint32_t            bitmap:30;
        uint32_t            lock  :1;
        uint32_t            alt   :1;
        uint16_t            min_idx;
        uint16_t            max_idx;
    } v;
} bleafMeta;

typedef struct leafEntry {
    float    k;
    uint32_t ch;
};

/**
 * bleaf: leaf node
 *
 * We guarantee that each leaf must have >=1 key.
 */
class bleaf {
public:
    uint32_t            bitmap:30;
    uint32_t            lock  :1;
    uint32_t            alt   :1;
    uint16_t            min_idx;
    uint16_t            max_idx;
    leafEntry           ent[LEAF_KEY_NUM];
    uint64_t            next[2];    // 0 next, 1 last

public:
    float  & k(int idx)  { return ent[idx].k; }
    uint32_t & ch(int idx) { return ent[idx].ch; }


    int num() {return countBit(bitmap);}
    bleaf * nextSibling() {
        return (bleaf*)nvm_allocator->getVirtualAddr(next[0]); 
    }

    void setBothHalfWord(bleafMeta *m) {
        bleafMeta * my_meta = (bleafMeta *)this;
        my_meta->halfWord4B[1] = m->halfWord4B[1];
        my_meta->halfWord4B[0] = m->halfWord4B[0];
    }

    void setlastLeaf(bleaf* vaddr) {
        next[1] = nvm_allocator->getPPointer(vaddr);
    }

    void setNextLeaf(bleaf* vaddr) {
        next[0] = nvm_allocator->getPPointer(vaddr);
    }

    bleaf* getLastLeaf() {
        return (bleaf*)nvm_allocator->getVirtualAddr(this->next[1]);
    }

    bleaf* getNextLeaf() {
        return (bleaf*)nvm_allocator->getVirtualAddr(this->next[0]);
    }

    float& getMinv() {
        return ent[min_idx].k;
    }

    float& getMaxv() {
        return ent[max_idx].k;
    }
}; // bleaf

/* ---------------------------------------------------------------------- */

class treeMeta {
 public:
    int        root_level; // leaf: level 0, parent of leaf: level 1
    Pointer8B  tree_root;
    uint64_t*   first_leaf; // on NVM
    uint64_t*  leaf_num;   // on NVM

 public:
    treeMeta(void *nvm_address, bool recover=false)
    { 
        root_level = 0; 
        tree_root=NULL; 
        first_leaf= (uint64_t*)nvm_address;
        leaf_num = (uint64_t*)(nvm_address + sizeof(uint64_t));
        if (! recover) {
            setFirstLeaf(NULL);
            setLeafNum(0);
        }
    }

    void setFirstLeaf(bleaf * leaf)
    {
         *first_leaf= nvm_allocator->getPPointer(leaf);
         clwb(first_leaf); sfence();
    }

    void setLeafNum(uint64_t num) {
        *leaf_num = num;
        clwb(leaf_num); sfence();
    }
 
}; // treeMeta


/* ---------------------------------------------------------------------- */

class lbtree: public tree {
  public:  // root and level
    
    treeMeta * tree_meta;
    
  public:
    lbtree(void *nvm_address, bool recover=false)
    {tree_meta= new treeMeta(nvm_address, recover);
     if (!tree_meta) {perror("new"); exit(1);}
    }

    ~lbtree()
    {delete tree_meta;}

  private:
    
    int bulkloadSubtree(KV *input, int start_key, int num_key, 
                        float bfill, int target_level,
                        Pointer8B pfirst[], int n_nodes[], uint64_t& g_memory);
    
    int bulkloadToptree(Pointer8B ptrs[], key_type keys[], int num_key,
                        float bfill, int cur_level, int target_level,
                        Pointer8B pfirst[], int n_nodes[]);

    void getMinMaxKey(bleaf *p, key_type &min_key, key_type &max_key);

    void getKeyPtrLevel(Pointer8B pnode, int pnode_level, key_type left_key,
         int target_level, Pointer8B ptrs[], key_type keys[], int &num_nodes,
         bool free_above_level_nodes);

    // sort pos[start] ... pos[end] (inclusively)
    void qsortBleaf(bleaf *p, int start, int end, int pos[]);

  public:
    int bulkload (int keynum, KV *input, float bfill, uint64_t& g_memory);
    
    // bulkload and recover the tree in nvm and return root level   
    int reBulkload(uint64_t& g_memory);

    // given a search key, perform the search operation
    // return the leaf node pointer and the position within leaf node
    void * lookup (key_type key, int *pos);
    
    void * get_recptr (void *p, int pos)
    {
        return (void*)((bleaf *)p)->ch(pos);
    }
    
private:
    void print (Pointer8B pnode, int level);
    void check (Pointer8B pnode, int level, key_type &start, key_type &end, bleaf* &ptr);
    void checkFirstLeaf(void);
    
public:
    void print ()
    {
        print (tree_meta->tree_root, tree_meta->root_level);
    }
    
    void check (key_type *start, key_type *end)
    {
        bleaf *ptr= NULL;
        check (tree_meta->tree_root, tree_meta->root_level, *start, *end, ptr);
        checkFirstLeaf();
    }
    
    int level () {return tree_meta->root_level;}
    
}; // lbtree

/* ---------------------------------------------------------------------- */
#endif /* _LBTREE_H */
