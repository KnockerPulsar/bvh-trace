// https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/
#include "vec_math.h"
#include "random.h"
#include "PPMImage.h"
#include "bvh_node.h"
#include "bvh.h"
#include "tri.h"
#include "ray.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <ostream>
#include <type_traits>
#include <chrono>

#define RES_X 640.0f
#define RES_Y 640.0f

extern Tri tri[N];
extern uint triIdx[N];

extern BVHNode bvhNode[N*2 - 1];
extern uint rootNodeIdx, nodesUsed; // Root node is used by default;
                                     
void ReadUnityModel() {
  FILE* file = fopen( "assets/unity.tri", "r" );
  float a, b, c, d, e, f, g, h, i;
  for (int t = 0; t < N; t++) 
  {
    if(fscanf( file, "%f %f %f %f %f %f %f %f %f\n", 
          &a, &b, &c, &d, &e, &f, &g, &h, &i ) > 0) {
      tri[t].vertex0 = float3( a, b, c );
      tri[t].vertex1 = float3( d, e, f );
      tri[t].vertex2 = float3( g, h, i );
    }
  }
  fclose( file );
}

void IntersectTri(Ray& ray, const Tri& tri) {
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

float IntersectAABB(  Ray& ray, const float3 bmin, const float3 bmax ) {
  float tx1 = (bmin.x - ray.O.x) / ray.D.x; float tx2 = (bmax.x - ray.O.x) / ray.D.x;

  float tmin = fmin(tx1, tx2), tmax = fmax(tx1, tx2);

  float ty1 = (bmin.y - ray.O.y) / ray.D.y;
  float ty2 = (bmax.y - ray.O.y) / ray.D.y;

  tmin = fmax(tmin, fmin(ty1, ty2)), tmax = fmin(tmax, fmax(ty1, ty2));

  float tz1 = (bmin.z - ray.O.z) / ray.D.z;
  float tz2 = (bmax.z - ray.O.z) / ray.D.z;

  tmin = fmax(tmin, fmin(tz1, tz2)), tmax = fmin(tmax, fmax(tz1, tz2));
  bool hit = tmax >= tmin && tmin < ray.t && tmax > 0;
  return hit? tmin : 1e30f;
}

void IntersectBVH(Ray& ray, const uint nodeIndex) {
  BVHNode* node = &bvhNode[nodeIndex];
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

    float dist1 = IntersectAABB(ray, child1->aabbMin, child1->aabbMax);
    float dist2 = IntersectAABB(ray, child2->aabbMin, child2->aabbMax);

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

int main() {

  Image output = Image::fromColor(RES_X, RES_Y, float3(0));

  ReadUnityModel();
  BuildBVH();

  float3 p0(-1, 1, 2), p1(1, 1, 2), p2(-1, -1, 2);
  Ray ray;
  ray.O = float3(-1.5f, -0.2f, -2.5f);

  auto start = std::chrono::high_resolution_clock::now();
  for ( int y = 0; y < RES_Y; y++) {
    // printf("%.2f\r", ((float)y+1)/RES_Y * 100.0f);
    // std::flush(std::cout);
    for ( int x = 0; x < RES_X; x++ ) {
      /*
                   |
             p0    |     p1
                   |
         ----------|------------
                   |     
             p2    |
                   |

         This enables us to get canvas positions in the square formed by
         p0, p1, p2.
     */

      float3 pixelPos = ray.O + p0 + (p1-p0) * (x/RES_X) + (p2-p0) * (y/RES_Y);
      ray.D = normalize(pixelPos - ray.O);
      ray.t = 1e30f;

#if 0
      for ( int i = 0; i < N; i++) 
        IntersectTri(ray, tri[i]);
#else 
      IntersectBVH(ray, rootNodeIdx);
#endif

      uint c = (255 - ray.t * 75);

      if(ray.t < 1e30f) {
        output.writeToPixel(x, y, float3(c));      
      }
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

  printf("Tracing took %ld milliseconds\n", duration); 
  output.save("./output.ppm");

  return 0;
}

