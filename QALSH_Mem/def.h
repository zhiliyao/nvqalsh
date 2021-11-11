#ifndef __DEF_H
#define __DEF_H

// -----------------------------------------------------------------------------
//  Macros
// -----------------------------------------------------------------------------
#define MIN(a, b)	(((a) < (b)) ? (a) : (b))
#define MAX(a, b)	(((a) > (b)) ? (a) : (b))
#define SQR(x)		((x) * (x))
#define SUM(x, y)	((x) + (y))
#define DIFF(x, y)	((y) - (x))
#define SWAP(x, y)	{ int tmp=x; x=y; y=tmp; }

// -----------------------------------------------------------------------------
//  Constants
// -----------------------------------------------------------------------------
const int   TOPK[]     = { 100 };
const int   MAX_ROUND  = 1;
const int   MAXK       = 100;

const int   SCAN_SIZE       = 128;
const int   SCAN_SIZE_TREE  = 128;
const int   GROUP_SIZE      = 64;
const int   CANDIDATES      = 100;

const float MAXREAL    = 3.402823466e+38F;
const float MINREAL    = -MAXREAL;
const int   MAXINT     = 2147483647;
const int   MININT     = -MAXINT;

const int   SIZEBOOL   = (int) sizeof(bool);
const int   SIZEINT    = (int) sizeof(int);
const int   SIZECHAR   = (int) sizeof(char);
const int   SIZEFLOAT  = (int) sizeof(float);
const int   SIZEDOUBLE = (int) sizeof(double);

const float E          = 2.7182818F;
const float PI         = 3.141592654F;
const float FLOATZERO  = 1e-6F;
const float ANGLE      = PI / 8.0f;

const uint64_t   POOL_SIZE      = 1900000000;
const uint64_t   DATA_POOL_SIZE = 1900000000;
const uint64_t   MEM_POOL_SIZE  = 1024 * 1024 * 1024;  // 1 GB

#endif // __DEF_H
