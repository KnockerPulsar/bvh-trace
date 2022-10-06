#include "defs.h"
#include "vec_math.h"
#include "tlas_node.h"
#include "bvh.h"
#include "bvh_instance.h"

class TLAS {
  public:
    TLAS() = default;
    TLAS(BVHInstance* bvhList, int n);
    void Build();
    void Intersect( bvt::Ray & ray );
  
  private:
    int FindBestMatch( int* list, int N, int A);

    TLASNode* tlasNode = nullptr;
    BVHInstance* blas = nullptr;
    uint nodesUsed, blasCount;
};
