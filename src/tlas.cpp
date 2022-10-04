#include "tlas.h"
#include "bvh.h"
#include "vec_math.h"
#include "util.h"

#include <array>
#include <cstdlib>
#include <fstream>

TLAS tlas;

TLAS::TLAS( bvt::BVH* bvhList, int n) {
  blas = bvhList;
  blasCount = n;

  tlasNode = (TLASNode*) aligned_alloc(64, sizeof(TLASNode) * 2 * n);
  nodesUsed = 2; // Skip one for cache alignment
}

void TLAS::Build() {

  tlasNode[2].leftBLAS = 0;
  tlasNode[2].aabbMin = make_float3(-100);
  tlasNode[2].aabbMax = make_float3(+100);
  tlasNode[2].isLeaf = true;

  tlasNode[3].leftBLAS = 1;
  tlasNode[3].aabbMin = make_float3(-100);
  tlasNode[3].aabbMax = make_float3(+100);
  tlasNode[3].isLeaf = true;

  tlasNode[0].leftBLAS = 2;
  tlasNode[0].aabbMin = make_float3(-100);
  tlasNode[0].aabbMax = make_float3(+100);
  tlasNode[0].isLeaf = false;
}

void TLAS::Intersect(bvt::Ray &ray) {
  TLASNode* node = &tlasNode[0], *stack[64];
  uint stackPtr = 0;

  while (1) {
    if(node->isLeaf) {
      blas[node->leftBLAS].Intersect(ray);
      if(stackPtr == 0) break;
      else              node = stack[--stackPtr];

      continue;
    }

    TLASNode* child1 = &tlasNode[node->leftBLAS];
    TLASNode* child2 = &tlasNode[node->leftBLAS + 1];

    float dist1 = INTERSECT_AABB(ray, child1);
    float dist2 = INTERSECT_AABB(ray, child2);

    if(dist1 > dist2) {
      std::swap(dist1, dist2);
      std::swap(child1, child2);
    }

    if(dist1 == 1e30f) {
      if(stackPtr == 0) break;
      else              node = stack[--stackPtr];
    } else {
      node = child1;
      if(dist2 != 1e30f) stack[stackPtr++] = child2;
    }
  }
}
