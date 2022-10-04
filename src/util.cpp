#include "util.h"
#include "bvh.h"


  template<unsigned i>
float vecElem( __m128 V)
{
  // shuffle V so that the element that you want is moved to the least-
  // significant element of the vector (V[0])
  V = _mm_shuffle_ps(V, V, _MM_SHUFFLE(i, i, i, i));
  // return the value in V[0]
  return _mm_cvtss_f32(V);
}

float IntersectAABB(  bvt::Ray& ray, const float3 bmin, const float3 bmax ) {
  float tx1 = (bmin.x - ray.O.x) * ray.rD.x; float tx2 = (bmax.x - ray.O.x) * ray.rD.x;

  float tmin = fmin(tx1, tx2), tmax = fmax(tx1, tx2);

  float ty1 = (bmin.y - ray.O.y) * ray.rD.y;
  float ty2 = (bmax.y - ray.O.y) * ray.rD.y;

  tmin = fmax(tmin, fmin(ty1, ty2)), tmax = fmin(tmax, fmax(ty1, ty2));

  float tz1 = (bmin.z - ray.O.z) * ray.rD.z;
  float tz2 = (bmax.z - ray.O.z) * ray.rD.z;

  tmin = fmax(tmin, fmin(tz1, tz2)), tmax = fmin(tmax, fmax(tz1, tz2));
  bool hit = tmax >= tmin && tmin < ray.t && tmax > 0;
  return hit? tmin : 1e30f;
}

float IntersectAABB_SSE(  bvt::Ray& ray, const __m128 bmin4, const __m128 bmax4 ) { 
  static __m128 mask4 = _mm_cmpeq_ps(_mm_setzero_ps(), _mm_set_ps(1,0,0,0));

  __m128 t1 = VMUL(VSUB(VAND(bmin4, mask4), ray.O4), ray.rD4);
  __m128 t2 = VMUL(VSUB(VAND(bmax4, mask4), ray.O4), ray.rD4);

  __m128 vmax4 = VMAX(t1,t2), vmin4 = VMIN(t1,t2);
  float tmax = std::min(vecElem<0>(vmax4), std::min(vecElem<1>(vmax4), vecElem<2>(vmax4)));
  float tmin = std::max(vecElem<0>(vmin4), std::max(vecElem<1>(vmin4), vecElem<2>(vmin4)));

  bool hit = (tmax >= tmin) && tmin < ray.t && tmax > 0;
  return hit? tmin : 1e30f;
}
