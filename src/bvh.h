#pragma once

#include "aabb.h"
#include "defs.h"
#include "bvh_node.h"
#include "ray.h"
#include "tri.h"

namespace bvt {
  class BVH {
    public:
      BVHNode* bvhNode = nullptr;

      BVH() = default;
      BVH(char* triFile, int n);
      void Build();
      void Refit();
      void Intersect(bvt::Ray& ray);
    private:
      void Subdivide(uint nodeIdx);
      float EvaluateSAH(BVHNode& node, int axis, float pos);
      void UpdateNodeBounds(uint nodeIdx);
      float FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos);
      float CalculateNodeCost(BVHNode& node);

      Tri* tri = nullptr;
      uint* triIdx = nullptr;
      uint nodesUsed = 2, triCount = 0;
  };
}

