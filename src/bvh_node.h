#pragma once

#include <defs.h>
#include <vec_math.h>
#include <xmmintrin.h>

struct BVHNode {
  union {
    struct {
      float3 aabbMin;
      // How we interpret this mysterious ‘leftFirst’ variable depends on
      // triCount. If it is 0, leftFirst contains the index of the left 
      // child node. Otherwise, it contains the index of the first triangle 
      // index.
      union { uint leftNode, firstTriIndex; };
    };
    __m128 aabbMin4;
  };
  union {
    struct { float3 aabbMax; uint triCount; };
    __m128 aabbMax4;
  };

  bool isLeaf() { return triCount > 0; }
};
