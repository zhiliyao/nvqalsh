/**
 *
 * @file tree.h
 * @author  Shimin Chen <shimin.chen@gmail.com>, Jihang Liu, Leying Chen
 * @version 1.0
 *
 * @section LICENSE
 *
 * TBD
 *
 * @section DESCRIPTION
 *
 * The tree class defines the methods of trees.
 */

#ifndef _BTREE_TREE_H
#define _BTREE_TREE_H
/* ---------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <thread>
#include <atomic>

// #include <immintrin.h>
/* ---------------------------------------------------------------------- */
/*                            Default Parameters                          */
/* ---------------------------------------------------------------------- */

// the size of a tree node
#define NONLEAF_LINE_NUM        4    // 256B
#define LEAF_LINE_NUM           4    // 256B

// the number of leaf nodes to prefetch ahead in jump pointer array 
// prefetching
#ifndef PREFETCH_NUM_AHEAD
#define PREFETCH_NUM_AHEAD	3
#endif

/* ---------------------------------------------------------------------- */
/*                 Node Size, Key Size, and Pointer Size                  */
/* ---------------------------------------------------------------------- */

// node size
#define NONLEAF_SIZE    (CACHE_LINE_SIZE * NONLEAF_LINE_NUM)
#define LEAF_SIZE       (CACHE_LINE_SIZE * LEAF_LINE_NUM)

// key size and pointer size: 8B
typedef double key_type;
#define KEY_SIZE             8   /* size of a key in tree node */
#define POINTER_SIZE         8   /* size of a pointer/value in node */
#define ITEM_SIZE            8   /* key size or pointer size */
#include <limits.h>
#include <cfloat>
#define MAX_KEY		DBL_MAX
#define MIN_KEY		DBL_MIN

// #include "keyinput.h"
#include "nodepref.h"
#include "mempool.h"
#include "nvm-common.h"
#include "util.h"
// #include "performance.h"

/* ---------------------------------------------------------------------- */
/*                            Useful funcions                             */
/* ---------------------------------------------------------------------- */
/* GCC builtin functions

int __builtin_ffs (unsigned int x)
Returns one plus the index of the least significant 1-bit of x, or if x is
zero, returns zero.

int __builtin_popcount (unsigned int x)
Returns the number of 1-bits in x.

*/

#define bitScan(x)  __builtin_ffs(x)
#define countBit(x) __builtin_popcount(x)

static inline unsigned long long rdtsc(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

#define Min(x,y)    ((x)<=(y) ? (x) : (y))
#define Max(x,y)    ((x)<=(y) ? (y) : (x))

// compute ceiling(x/y) and floor(x/y)
#define ceiling(x, y)  (((x) + (y) - 1) / (y))
#define floor(x, y)    ((x) / (y))

#define Swap(x, y) \
do { auto _t=(x); (x)=(y); (y)=_t; } while(0)

/* ---------------------------------------------------------------------- */
class tree {
 public:
   virtual int bulkload (int keynum, KV *input, float bfill, uint64_t& g_memory)
   {
	fprintf (stderr, "Not implemented!\n");
	exit (1);
	return 0;
   }

   virtual int reBulkload(uint64_t& g_memory) {}

  /**
   * randomize the key orders in nodes (only for unsorted/bitmap trees)
   */
   virtual void randomize() {}

  /**
   * given a search key, perform the search operation
   *
   * @param key   the search key
   * @param pos   the position to return
   * @return the leaf node to return
   *
   * If a match to the given search key is found, then the leaf node and the
   * matching key position is returned.  If a match is not found, then the leaf
   * node and the position to a previous key is returned.
   */
   virtual void * lookup (key_type key, int *pos)
   {
	fprintf (stderr, "Not implemented!\n");
	exit (1);
	return NULL;
   }

  /**
   * obtain the record pointer
   *
   * @param p     leaf node pointer
   * @param pos   index entry position in the leaf node
   * @return      the associated record pointer
   */
   virtual void * get_recptr (void *p, int pos)
   {
        fprintf (stderr, "Not implemented!\n");
	exit (1);
	return NULL;
   }

  /**
   * insert an index entry
   *
   * @param key   the index key
   * @param ptr   the record pointer
   */
   virtual void insert (key_type key, void * ptr)
   {
       fprintf (stderr, "Not implemented!\n");
       exit (1);
   }

  /**
   * delete an index entry
   *
   * @param key   the index key to delete
   */
   virtual void del (key_type key)
   {
	fprintf (stderr, "Not implemented!\n");
	exit (1);
   }

  /**
   * print the tree structure
   */
   virtual void print ()
   {
	fprintf (stderr, "Not implemented!\n");
	exit (1);
   }

  /**
   * check the correctness of the tree structure
   *
   * @param start   the start key of the tree
   * @param end     the end key of the tree
   */
   virtual void check (key_type *start, key_type *end)
   {
	fprintf (stderr, "Not implemented!\n");
	exit (1);
   }

  /**
   * obtain the level parameter of the tree
   */
   virtual int level ()
   {
        fprintf (stderr, "Not implemented!\n");
        exit (1);
	return 0;
   }

}; // tree

/* ---------------------------------------------------------------------- */
extern tree * the_treep;
extern int    worker_thread_num;
extern const char * nvm_file_name;

extern int parse_command (int argc, char **argv);

#ifdef INSTRUMENT_INSERTION
extern int insert_total;       // insert_total=
extern int insert_no_split;    //              insert_no_split
extern int insert_leaf_split;  //             +insert_leaf_split
extern int insert_nonleaf_split;//            +insert_nonleaf_split
extern int total_node_splits;  // an insertion may cause multiple node splits
#endif // INSTRUMENT_INSERTION

// a specific tree should implement a child class of tree
// and the following function
extern tree * initTree(void *nvm_addr, bool recover);

/* ---------------------------------------------------------------------- */
#endif /* _BTREE_TREE_H */
