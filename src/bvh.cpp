#include "bvh.h"
#include "aabb.h"
#include "bvh_node.h"
#include <numeric>
#include <stdio.h>

#define EVALUATE_SAH 0

Tri tri[N];
uint triIdx[N];

BVHNode bvhNode[N*2 - 1];
uint rootNodeIdx = 0, nodesUsed = 1; // Root node is used by default;

void UpdateNodeBounds(uint nodeIndex) {
  BVHNode& node = bvhNode[nodeIndex];
  node.aabbMin = float3(1e30f);
  node.aabbMax = float3(-1e30f);

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

void Subdivide(uint nodeIndex) {
  BVHNode& node = bvhNode[nodeIndex];

  float3 extent = node.aabbMax - node.aabbMin;

  int axis = 0;
  float splitPos;

#if EVALUATE_SAH
  int bestAxis = -1;
  float bestPos = 0, bestCost = 1e30f;

  float parentArea = extent.x * extent.y + extent.x * extent.z + extent.y * extent.z;
  float parentCost = node.triCount * parentArea;

  for( int axis = 0; axis < 3; axis++) {
    for( uint i = 0; i < node.triCount; i++) {
      Tri& triangle = tri[triIdx[node.firstTriIndex + i]];
      float candidatePos = triangle.centroid[axis];
      float cost = EvaluateSAH(node, axis, candidatePos);
      if(cost < bestCost)
        bestPos = candidatePos, bestAxis = axis, bestCost = cost;
    }
  }

  axis = bestAxis;
  splitPos = bestPos;

  if(bestCost >= parentCost) return; 
#else
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
  printf("Building BVH\n");
  for (int i = 0; i < N; i++) {
    tri[i].centroid = (tri[i].vertex0 + tri[i].vertex1 + tri[i].vertex2) * 0.333f;
    triIdx[i]  = i;
  }

  BVHNode& root = bvhNode[rootNodeIdx];
  root.leftNode = 0;
  root.firstTriIndex = 0; 
  root.triCount = N;

  UpdateNodeBounds(rootNodeIdx);
  Subdivide(rootNodeIdx);
  printf("Finished building BVH\n");
}

