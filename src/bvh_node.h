#pragma once

#include <defs.h>
#include <vec_math.h>

struct BVHNode {
  float3 aabbMin, aabbMax;

  // How we interpret this mysterious ‘leftFirst’ variable depends on
  // triCount. If it is 0, leftFirst contains the index of the left 
  // child node. Otherwise, it contains the index of the first triangle 
  // index.
  union { uint leftNode, firstTriIndex; };
  uint triCount;

  bool isLeaf() { return triCount > 0; }
};
