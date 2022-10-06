#include "tlas.h"
#include "bvh.h"
#include "bvh_instance.h"
#include "vec_math.h"
#include "util.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>

TLAS tlas;

TLAS::TLAS( BVHInstance* bvhList, int n) {
  blas = bvhList;
  blasCount = n;

  tlasNode = (TLASNode*) aligned_alloc(64, sizeof(TLASNode) * 2 * n);
  nodesUsed = 2; // Skip one for cache alignment
}

int TLAS::FindBestMatch( int* list, int N, int A) {
  float smallest = 1e30f;
  int bestB = -1;

  for (int B = 0 ; B < N; B++) if (B != A) {
    float3 bmax = fmaxf(tlasNode[list[A]].aabbMax, tlasNode[list[B]].aabbMax);
    float3 bmin = fminf(tlasNode[list[A]].aabbMin, tlasNode[list[B]].aabbMin);

    float3 e = bmax - bmin;
    float surfaceArea = e.x * e.y + e.y * e.z + e.z * e.x;
    if(surfaceArea < smallest) smallest = surfaceArea, bestB = B;
  }

  return bestB;
}

void TLAS::Build() {
  int nodeIdx[256], nodeIndices = blasCount;
  nodesUsed = 1;

  for (uint i = 0; i < blasCount; i++) {
    nodeIdx[i] = nodesUsed;
    tlasNode[nodesUsed].aabbMin = blas[i].bounds.min;
    tlasNode[nodesUsed].aabbMax = blas[i].bounds.max;
    tlasNode[nodesUsed].BLAS = i;

    // Leaf node
    tlasNode[nodesUsed++].leftRight = 0;
  }

  int A = 0, B = FindBestMatch(nodeIdx, nodeIndices, A);
  while (nodeIndices > 1) {
    int C = FindBestMatch(nodeIdx, nodeIndices, B);
    if(A == C) {
      int nodeIdxA = nodeIdx[A], nodeIdxB = nodeIdx[B];

      TLASNode& nodeA = tlasNode[nodeIdxA];
      TLASNode& nodeB = tlasNode[nodeIdxB];

      TLASNode& newNode = tlasNode[nodesUsed];

      newNode.leftRight = nodeIdxA + (nodeIdxB << 16);
      newNode.aabbMin = fminf(nodeA.aabbMin, nodeB.aabbMin);
      newNode.aabbMax = fmaxf(nodeA.aabbMax, nodeB.aabbMax);

      nodeIdx[A] = nodesUsed++;
      nodeIdx[B] = nodeIdx[nodeIndices - 1];

      B = FindBestMatch(nodeIdx, --nodeIndices, A);
    } else {
      A = B, B = C;
    }

    tlasNode[0] = tlasNode[nodeIdx[A]];
  }
} 

void TLAS::Intersect(bvt::Ray &ray) {
  ray.rD = make_float3(1/ray.D.x, 1/ray.D.y, 1/ray.D.z);

  TLASNode* node = &tlasNode[0], *stack[64];
  uint stackPtr = 0;

  while (1) {
    if(node->isLeaf()) {
      blas[node->BLAS].Intersect(ray);
      if(stackPtr == 0) break;
      else              node = stack[--stackPtr];

      continue;
    }

    TLASNode* child1 = &tlasNode[node->leftRight & 0xffff];
    TLASNode* child2 = &tlasNode[node->leftRight >> 16];

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
