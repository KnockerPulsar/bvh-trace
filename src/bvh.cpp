#include "bvh.h"
#include "aabb.h"
#include "bvh_node.h"
#include <chrono>
#include <numeric>
#include <stdio.h>

#define EVALUATE_SAH 1
#define SAH_TEST_PLANES 4

Tri tri[N];
uint triIdx[N];

BVHNode bvhNode[N*2 - 1];
uint rootNodeIdx = 0, nodesUsed = 1; // Root node is used by default;

void UpdateNodeBounds(uint nodeIndex) {
  BVHNode& node = bvhNode[nodeIndex];
  node.aabbMin = make_float3(1e30f);
  node.aabbMax = make_float3(-1e30f);

  for (uint first = node.firstTriIndex, i = 0; i < node.triCount; i++) {
    uint leafTriIndex = triIdx[first + i];
    Tri& leafTri = tri[leafTriIndex];

    node.aabbMin = fminf(node.aabbMin, fminf(leafTri.vertex0, fminf(leafTri.vertex1, leafTri.vertex2)));
    node.aabbMax = fmaxf(node.aabbMax, fmaxf(leafTri.vertex0, fmaxf(leafTri.vertex1, leafTri.vertex2)));
    // node.aabbMax = node.aabbMax * 1.13f; // For some reason the BVH don't fit quite right?
  }
}

float EvaluateSAH( BVHNode& node, int axis, float pos) {
  AABB leftBox, rightBox; 
  int leftCount = 0, rightCount = 0;

  for (uint i = 0; i < node.triCount; i++) {
    Tri& triangle = tri[triIdx[node.firstTriIndex + i]];
    if(triangle.centroid[axis] < pos) {
      leftCount ++;
      leftBox.grow(triangle.vertex0);
      leftBox.grow(triangle.vertex1);
      leftBox.grow(triangle.vertex2);
    } else {
      rightCount ++;
      rightBox.grow(triangle.vertex0);
      rightBox.grow(triangle.vertex1);
      rightBox.grow(triangle.vertex2);
    }
  }
  float cost = leftCount * leftBox.area() + rightCount * rightBox.area();
  return cost > 0? cost: 1e30f;
}

float FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos) {
  float bestCost = 1e30f;

  for( int a = 0; a < 3; a++) {
    float boundsMin = node.aabbMin[a];
    float boundsMax = node.aabbMax[a];

    if(boundsMin == boundsMax) continue;
    float scale = (boundsMax - boundsMin) / SAH_TEST_PLANES;

    for( uint i = 1; i < SAH_TEST_PLANES; i++) {
      float candidatePos = boundsMin + i * scale;
      float cost = EvaluateSAH(node, a, candidatePos);
      if(cost < bestCost)
        splitPos = candidatePos, axis = a, bestCost = cost;
    }
  }

  return bestCost;
}

float CalculateNodeCost( BVHNode& node ) {
  float3 e = node.aabbMax - node.aabbMin;
  float surfaceArea = e.x * e.y + e.x * e.z + e.y * e.z;
  return node.triCount * surfaceArea;
}

void Subdivide(uint nodeIndex) {
  BVHNode& node = bvhNode[nodeIndex];

  int axis = 0;
  float splitPos;

#if EVALUATE_SAH
  float bestCost = FindBestSplitPlane(node, axis, splitPos);
  float noSplitCost = CalculateNodeCost(node);
  if(bestCost >= noSplitCost) return; 
#else
  float3 extent = node.aabbMax - node.aabbMin;
  if(extent.y > extent.x) axis = 1;
  if(extent.z > extent[axis]) axis = 2;

  splitPos = node.aabbMin[axis] + extent[axis] * 0.5f;
#endif

  int i = node.leftNode;
  int j = i + node.triCount - 1;

  while (i <= j) {
    if (tri[triIdx[i]].centroid[axis] < splitPos) {
      i++;
    } else {
      std::swap(triIdx[i], triIdx[j--]);
    }
  }

  int leftCount = i - node.firstTriIndex;
  if(leftCount == 0 || leftCount == node.triCount) return;

  int leftChildIndex = nodesUsed++;
  int rightChildIndex = nodesUsed++;


  bvhNode[leftChildIndex].firstTriIndex = node.firstTriIndex;
  bvhNode[leftChildIndex].triCount = leftCount;

  bvhNode[rightChildIndex].firstTriIndex = i;
  bvhNode[rightChildIndex].triCount = node.triCount - leftCount;

  node.leftNode = leftChildIndex;
  node.triCount = 0;

  UpdateNodeBounds(leftChildIndex);
  UpdateNodeBounds(rightChildIndex);

  Subdivide(leftChildIndex);
  Subdivide(rightChildIndex);
}

void BuildBVH() {
  for (int i = 0; i < N; i++) {
    tri[i].centroid = (tri[i].vertex0 + tri[i].vertex1 + tri[i].vertex2) * 0.333f;
    triIdx[i]  = i;
  }

  BVHNode& root = bvhNode[rootNodeIdx];
  root.leftNode = 0;
  root.firstTriIndex = 0; 
  root.triCount = N;

  UpdateNodeBounds(rootNodeIdx);

  printf("Building BVH\n");
  auto start = std::chrono::high_resolution_clock::now();
  Subdivide(rootNodeIdx);
  auto end = std::chrono::high_resolution_clock::now();
  long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
  printf("Finished building BVH, took %ld milliseconds\n", duration);
}

