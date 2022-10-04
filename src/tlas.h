#include "defs.h"
#include "vec_math.h"
#include "tlas_node.h"
#include "bvh.h"

class TLAS {
  public:
    TLAS() = default;
    TLAS( bvt::BVH* bvhList, int n);
    void Build();
    void Intersect( bvt::Ray & ray );
  
  private:
    TLASNode* tlasNode = nullptr;
    bvt::BVH* blas = nullptr;
    uint nodesUsed, blasCount;
};
