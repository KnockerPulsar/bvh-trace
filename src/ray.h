#pragma once
#include "vec_math.h"
#include <xmmintrin.h>

namespace bvt {
  struct Ray {
    Ray() { O4 = D4 = rD4 = _mm_set1_ps(1); };
    union { struct { float3 O; float dummy1; }; __m128 O4; };
    union { struct { float3 D; float dummy2; }; __m128 D4; };
    union { struct { float3 rD; float dummy3; }; __m128 rD4; };
    float t = 1e30f; 
  };
}
