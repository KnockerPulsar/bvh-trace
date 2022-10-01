#include "bvh.h"

Tri tri[N];
uint triIdx[N];

BVHNode bvhNode[N*2 - 1];
uint rootNodeIdx = 0, nodesUsed = 1; // Root node is used by default;
                                     //
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

void Subdivide(uint nodeIndex) {
  BVHNode& node = bvhNode[nodeIndex];
  if(node.triCount <= 2) return;

  float3 extent = node.aabbMax - node.aabbMin;

  int axis = 0;
  if(extent.y > extent.x) axis = 1;
  if(extent.z > extent[axis]) axis = 2;

  float splitPos = node.aabbMin[axis] + extent[axis] * 0.5f;

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
  Subdivide(rootNodeIdx);

}

