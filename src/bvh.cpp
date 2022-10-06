#include "bvh.h"
#include "aabb.h"
#include "bin.h"
#include "bvh_instance.h"
#include "bvh_node.h"
#include "tri.h"
#include "vec_math.h"
#include "util.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <numeric>
#include <raylib.h>
#include <stdio.h>

BVHInstance bvhInstance[256];

namespace bvt {
  void IntersectTri(bvt::Ray& ray, const Tri& tri) {
    float3 edge1 = tri.vertex1 - tri.vertex0;
    float3 edge2 = tri.vertex2 - tri.vertex0;

    const float3 h = cross(ray.D, edge2);
    const float det = dot(edge1, h);

    // Ray parallel to triangle
    if(det > -0.0001f  && det < 0.0001f) return;

    const float invDet = 1.0f/det;
    const float3 v0RayOrigin = ray.O - tri.vertex0;

    // First barycentric coordinate
    const float u = invDet * dot(v0RayOrigin, h);

    // Cant be inside the triangle if u < 0 or > 1
    if( u < 0 || u > 1.0f) return;

    const float3 q = cross(v0RayOrigin, edge1);

    // Second barycentric coordinate
    const float v = invDet * dot(ray.D, q);

    // Same for v
    // Notice that we check if u+v > 1.0f
    // This is because at the point of checking both u > 0 and v > 0
    // But we need to also check if their sum isn't > 1.0f
    if( v < 0 || u + v > 1.0f) return;

    // At this point we're certain we hit the triangle.
    //
    // How far along the ray we intersected the triangle.
    const float t = invDet * dot(edge2, q);

    // Account for numerical precision?
    if( t > 0.0001f ) {
      // Update the ray's intersection distance if a closer one is found.
      ray.t = std::min(ray.t, t);
    }

  }

  BVH::BVH(char * triFile, int n) {
    FILE* file = fopen( triFile, "r" );
    float a, b, c, d, e, f, g, h, i;

    tri = new Tri[n];
    triIdx = new uint[n];

    bvhNode = (BVHNode* )aligned_alloc(64, sizeof(BVHNode) * 2 * n);
    bvhNode[0].triCount = triCount = n;
    bvhNode[0].firstTriIndex = 2;


    float3 center = make_float3(0);
    for (int t = 0; t < n; t++) 
    {
      if(fscanf( file, "%f %f %f %f %f %f %f %f %f\n", 
            &a, &b, &c, &d, &e, &f, &g, &h, &i ) > 0) {
        tri[t].vertex0 = make_float3( a, b, c );
        tri[t].vertex1 = make_float3( d, e, f );
        tri[t].vertex2 = make_float3( g, h, i );

        center.x += a+d+g;
        center.y += b+e+h;
        center.z += c+f+i;
      }
    }


    center = center * (1.0f/n);

    /* printf("Center at (%f, %f, %f)\n", center.x, center.y, center.z); */

    fclose( file );

    Build();
  }

  void BVH::UpdateNodeBounds(uint nodeIndex) {
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

  float BVH::EvaluateSAH( BVHNode& node, int axis, float pos) {
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

  float BVH::FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos) {
    float bestCost = 1e30f;
    for( int a = 0; a < 3; a++) {

      float boundsMin = 1e30f, boundsMax = -1e30f; 
      for (int i = 0; i < node.triCount; i++) {
        Tri& triangle = tri[triIdx[ node.firstTriIndex + i ]];
        boundsMin = fminf(boundsMin, triangle.centroid[a]);
        boundsMax = fmaxf(boundsMax, triangle.centroid[a]);
      }

      if(boundsMin == boundsMax) continue;

      Bin bin[BINS];
      float scale = BINS / (boundsMax - boundsMin);

      for (uint i = 0; i < node.triCount; i++) {
        Tri& triangle = tri[triIdx[node.firstTriIndex + i]];

        int binIdx = std::min(BINS-1, (int)((triangle.centroid[a] - boundsMin) * scale)); 
        bin[binIdx].triCount++;
        bin[binIdx].bounds.grow(triangle.vertex0);
        bin[binIdx].bounds.grow(triangle.vertex1);
        bin[binIdx].bounds.grow(triangle.vertex2);
      }

      float leftArea[BINS-1], rightArea[BINS-1];
      int leftCount[BINS-1], rightCount[BINS-1];

      AABB leftBox, rightBox;
      int leftSum = 0, rightSum = 0;
      for (int i = 0; i < BINS - 1; i++) {
        leftSum += bin[i].triCount;
        leftCount[i] = leftSum;
        leftBox.grow(bin[i].bounds);
        leftArea[i] = leftBox.area();

        rightSum += bin[BINS - 1 - i].triCount;
        rightCount[BINS - 2 -i] = rightSum;
        rightBox.grow(bin[BINS - 1 - i].bounds);
        rightArea[BINS - 2 - i] = rightBox.area();
      }

      scale = (boundsMax - boundsMin) / BINS;
      for( uint i = 1; i < BINS - 1; i++) {
        float planeCost = leftCount[i] * leftArea[i] + rightCount[i] * rightArea[i];
        if(planeCost < bestCost)
          splitPos = boundsMin + scale * (i+1), axis = a, bestCost = planeCost;
      }
    }

    return bestCost;
  }

  float BVH::CalculateNodeCost( BVHNode& node ) {
    float3 e = node.aabbMax - node.aabbMin;
    float surfaceArea = e.x * e.y + e.x * e.z + e.y * e.z;
    return node.triCount * surfaceArea;
  }

  void BVH::Subdivide(uint nodeIndex) {
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

  void BVH::Build() {
    for (int i = 0; i < triCount; i++) {
      tri[i].centroid = (tri[i].vertex0 + tri[i].vertex1 + tri[i].vertex2) * 0.333f;
      triIdx[i]  = i;
    }

    for (int i = 1; i < triCount; i++) {
      bvhNode[i].triCount = 0;
      bvhNode[i].firstTriIndex = 0;
      bvhNode[i].aabbMin = make_float3(1e30f);
      bvhNode[i].aabbMax = make_float3(-1e30f);
    }

    UpdateNodeBounds(0);

    /* printf("Building BVH\n"); */
    /* auto start = std::chrono::high_resolution_clock::now(); */
    Subdivide(0);
    /* auto end = std::chrono::high_resolution_clock::now(); */
    /* long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count(); */
    /* printf("Finished building BVH, took %ld milliseconds\n", duration); */
  }

  // Bottom up refit
  // Should be used only for slight movements where triangles don't move between nodes
  void BVH::Refit() {
    for (int i = nodesUsed - 1; i > 1; i--) {
      BVHNode& node = bvhNode[i];

      if(node.isLeaf()) {
        UpdateNodeBounds(i);
        continue;
      }

      BVHNode& leftChild = bvhNode[node.leftNode];
      BVHNode& rightChild = bvhNode[node.leftNode + 1];

      node.aabbMin = fminf(leftChild.aabbMin, rightChild.aabbMin);
      node.aabbMax = fmaxf(leftChild.aabbMax, rightChild.aabbMax);
    }
  }

  void BVH::Intersect(bvt::Ray& ray) {
    BVHNode* node = &bvhNode[0];
    BVHNode* stack[64];
    uint stackPtr = 0;

    while(1) {
      if(node->isLeaf()) {
        for (uint i = 0; i < node->triCount; i++) {
          IntersectTri(ray, tri[triIdx[node->firstTriIndex + i]]);
        }

        if(stackPtr == 0) break;
        else              node = stack[--stackPtr];

        continue;
      }

      BVHNode* child1 = &bvhNode[node->leftNode];
      BVHNode* child2 = &bvhNode[node->leftNode + 1];

      float dist1 = INTERSECT_AABB(ray, child1); 
      float dist2 =  INTERSECT_AABB(ray, child2);  

      if(dist1 > dist2) {std::swap(dist1, dist2); std::swap(child1, child2); }
      if(dist1 == 1e30f) {
        if(stackPtr == 0) break;
        else              node = stack[--stackPtr];
      } else {
        node = child1;
        if(dist2 != 1e30f) stack[stackPtr++] = child2;
      }
    }
  }
}
