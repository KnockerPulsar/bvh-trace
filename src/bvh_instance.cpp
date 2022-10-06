#include "bvh_instance.h"
#include "vec_math.h"

void BVHInstance::Intersect(bvt::Ray &ray) {
  bvt::Ray backupRay = ray;

  ray.O = TransformPosition(ray.O, invTransform);
  ray.D = TransformVector(ray.D, invTransform);

  ray.rD = make_float3(1/ray.D.x, 1/ray.D.y, 1/ray.D.z);
  
  bvh->Intersect(ray);

  backupRay.t = ray.t;
  ray = backupRay;
}

void BVHInstance::SetTransform(const mat4& transform) {
  invTransform = transform.Inverted();
  float3 bmin = bvh->bvhNode[0].aabbMin, bmax = bvh->bvhNode[0].aabbMax;

  bounds = AABB();
  for (int i = 0; i < 8; i++) {
    bounds.grow(TransformPosition(make_float3(i & 1? bmax.x : bmin.x, i & 2? bmax.y : bmin.y, i & 4? bmax.z : bmin.z), transform));
  }

  // printf("Bounds min: (%f, %f, %f), Bounds max: (%f, %f, %f)\n", 
  //     bounds.min.x, bounds.min.y, bounds.min.z,
  //     bounds.max.x, bounds.max.y, bounds.max.z
  // );
}
