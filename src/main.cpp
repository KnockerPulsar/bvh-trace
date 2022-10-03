// https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/
#include "vec_math.h"
#include "random.h"
#include "PPMImage.h"
#include "bvh_node.h"
#include "bvh.h"
#include "tri.h"
#include "ray.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <ostream>
#include <type_traits>
#include <chrono>
#include <xmmintrin.h>
#include <raylib.h>
#include <rlgl.h>

#define RES_X 640.0f
#define RES_Y 640.0f
#define USE_SSE 1

#define VMUL(a,b) _mm_mul_ps(a,b)
#define VSUB(a,b) _mm_sub_ps(a,b)
#define VAND(a,b) _mm_and_ps(a,b)
#define VMAX(a,b) _mm_max_ps(a,b)
#define VMIN(a,b) _mm_min_ps(a,b)

  template<unsigned i>
float vecElem( __m128 V)
{
  // shuffle V so that the element that you want is moved to the least-
  // significant element of the vector (V[0])
  V = _mm_shuffle_ps(V, V, _MM_SHUFFLE(i, i, i, i));
  // return the value in V[0]
  return _mm_cvtss_f32(V);
}

extern Tri tri[N];
extern uint triIdx[N];

extern BVHNode bvhNode[N*2 - 1];
extern uint rootNodeIdx, nodesUsed; // Root node is used by default;

void ReadUnityModel() {
  FILE* file = fopen( "assets/bigben.tri", "r" );
  float a, b, c, d, e, f, g, h, i;

  float3 center = make_float3(0);
  for (int t = 0; t < N; t++) 
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

  center = center * (1.0f/N);

  printf("Center at (%f, %f, %f)\n", center.x, center.y, center.z);

  fclose( file );
}

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

void IntersectBVH(bvt::Ray& ray) {
  BVHNode* node = &bvhNode[rootNodeIdx];
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

#if USE_SSE
    float dist1 = IntersectAABB_SSE(ray, child1->aabbMin4, child1->aabbMax4);
    float dist2 = IntersectAABB_SSE(ray, child2->aabbMin4, child2->aabbMax4);
#else
    float dist1 = IntersectAABB(ray, child1->aabbMin, child1->aabbMax);
    float dist2 = IntersectAABB(ray, child2->aabbMin, child2->aabbMax);
#endif

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

void animate() {

}

int main() {

  bvt::Image output = bvt::Image::fromColor(RES_X, RES_Y, make_float3(0));

  ReadUnityModel();

  InitWindow(RES_X, RES_Y, "BVH Tracer");

  Texture frameBuffer = LoadTextureFromImage({
      .data = (void*)output.getData().data(),             
      .width = (int) RES_X,              
      .height = (int) RES_Y,
      .mipmaps = 1,
      .format = RL_PIXELFORMAT_UNCOMPRESSED_R8G8B8
      });

  BuildBVH();

  float3 p0{-1, 1, 2}, p1{1, 1, 2}, p2{-1, -1, 2};

  const int tile_size = 8;

  while(!WindowShouldClose()) {

    BuildBVH();

#pragma parallel for schedule(dynamic)
    auto start = std::chrono::high_resolution_clock::now();
    for (int tile = 0; tile < 6400; tile++) {
      int x = tile % 80, y = tile / 80;
      bvt::Ray ray;
      ray.O = make_float3(0, 3.5f, -4.5f);
      for (int v = 0; v < tile_size ; v++ ) for (int u = 0; u < tile_size; u++) {

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

        float3 pixelPos = ray.O + p0 + (p1-p0) * ((x * 8 +u)/RES_X) + (p2-p0) * ((y * 8+v)/RES_Y);
        ray.D = normalize(pixelPos - ray.O);
        ray.rD = make_float3(1/ray.D.x , 1/ray.D.y, 1/ray.D.z); 
        ray.t = 1e30f;

        IntersectBVH(ray);

        uint c = (255 - (int)((ray.t - 4) * 180));

        if(ray.t < 1e30f) {
          output.writeToPixel(x * 8 + u, y * 8 + v, make_float3(c));      
        }
      }
    }

    auto end = std::chrono::high_resolution_clock::now();
    long duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    printf("Tracing took %ld milliseconds\n", duration); 

    BeginDrawing();
    ClearBackground(BLACK);

    UpdateTexture(frameBuffer, (void*)output.getData().data());

    DrawTexture(frameBuffer, 0, 0, WHITE);
    EndDrawing();

    output.clear();
  }
  output.save("./output.ppm");

  return 0;
}

