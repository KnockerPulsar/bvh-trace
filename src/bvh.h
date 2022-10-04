#pragma once

#include "aabb.h"
#include "defs.h"
#include "bvh_node.h"
#include "ray.h"
#include "tri.h"

namespace bvt {
  class BVH {
    public:
      BVH() = default;
      BVH(char* triFile, int n);
      void Build();
      void Refit();
      void Intersect(bvt::Ray& ray);
      void SetTranform(const mat4& transform);
    private:
      void Subdivide(uint nodeIdx);
      float EvaluateSAH(BVHNode& node, int axis, float pos);
      void UpdateNodeBounds(uint nodeIdx);
      float FindBestSplitPlane(BVHNode& node, int& axis, float& splitPos);
      float CalculateNodeCost(BVHNode& node);

      BVHNode* bvhNode = nullptr;
      Tri* tri = nullptr;
      uint* triIdx = nullptr;
      uint nodesUsed = 2, triCount = 0;

      mat4 invTransform;
      AABB bounds;
  };
}

