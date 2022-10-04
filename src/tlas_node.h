#include "defs.h"
#include "vec_math.h"

#include <xmmintrin.h>

struct TLASNode {
  union { 
    struct {float3 aabbMin; uint leftBLAS;};
    __m128 aabbMin4;
  };
  union {
    struct {float3 aabbMax; uint isLeaf;};
    __m128 aabbMax4;
  };
};
