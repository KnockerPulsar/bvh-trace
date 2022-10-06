#include "config.h"
#include "bvh.h"
#include <cstdio>
#include <memory.h>

#define VMUL(a,b) _mm_mul_ps(a,b)
#define VSUB(a,b) _mm_sub_ps(a,b)
#define VAND(a,b) _mm_and_ps(a,b)
#define VMAX(a,b) _mm_max_ps(a,b)
#define VMIN(a,b) _mm_min_ps(a,b)

float IntersectAABB(  bvt::Ray& ray, const float3 bmin, const float3 bmax );
float IntersectAABB_SSE(  bvt::Ray& ray, const __m128 bmin4, const __m128 bmax4 );

#if USE_SSE 
#define INTERSECT_AABB(ray, node) IntersectAABB_SSE(ray, node->aabbMin4, node->aabbMax4) 
#else 
#define INTERSECT_AABB(ray, node) IntersectAABB(ray, node->aabbMin, node->aabbMax) 
#endif
