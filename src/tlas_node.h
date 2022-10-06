#include "defs.h"
#include "vec_math.h"

#include <xmmintrin.h>

struct TLASNode {
  union { 
    struct {
      float3 aabbMin; 

      // Left = first 16 bits
      // Right = second 16 bits
      uint leftRight;
    };
    __m128 aabbMin4;
  };
  union {
    struct {float3 aabbMax; uint BLAS;};
    __m128 aabbMax4;
  };

  bool isLeaf() { return leftRight == 0; }
};
