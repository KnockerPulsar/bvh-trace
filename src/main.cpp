// https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/
#include "vec_math.h"
#include "random.h"
#include "PPMImage.h"
#include <cmath>
#include <math.h>

#define N 12

struct Tri { float3 vertex0, vertex1, vertex2,  centroid ;};
struct Ray { float3 O, D; float t = 1e30f; };

Tri tri[N];

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
  const float u = dot(v0RayOrigin, h);

  // Cant be inside the triangle if u < 0 or > 1
  if( u < 0 || u > 1.0f) return;

  const float3 q = cross(v0RayOrigin, edge1);

  // Second barycentric coordinate
  const float v = dot(ray.D, q);

  // Same for v
  // Notice that we check if u+v > 1.0f
  // This is because at the point of checking both u > 0 and v > 0
  // But we need to also check if their sum isn't > 1.0f
  if( v < 0 || u + v > 1.0f) return;
  
  // At this point we're certain we hit the triangle.
  //
  // How far along the ray we intersected the triangle.
  const float t = dot(edge2, q);

  // Account for numerical precision?
  if( t > 0.0001f )
    // Update the ray's intersection distance if a closer one is found.
    ray.t = std::min(ray.t, t);

}


int main() {

  Image output = Image::magenta(640, 640);

  for (int i = 0 ; i < N; i++) {
    float3 r0 (RandomFloat(), RandomFloat(), RandomFloat()); 
    float3 r1 (RandomFloat(), RandomFloat(), RandomFloat()); 
    float3 r2 (RandomFloat(), RandomFloat(), RandomFloat()); 

    tri[i].vertex0 = r0 * 9 - float3(5);
    tri[i].vertex1 = tri[i].vertex0 + r1;
    tri[i].vertex2 = tri[i].vertex0 + r2;
  }

  float3 camPos(0, 0, -18);
  float3 p0(-1, 1, -15), p1(1, 1, -15), p2(-1, -1, -15);
  Ray ray;

        
  for ( int y = 0; y < 640; y++) {
    for ( int x = 0; x < 640; x++ ) {
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
      float3 pixelPos = p0 + (p1-p0) * (x/640.0f) + (p2-p0) * (y/640.0f);
      ray.O = camPos;
      ray.D = normalize(pixelPos - ray.O);
      ray.t = 1e30f;

      for ( int i = 0; i < N; i++) 
        IntersectTri(ray, tri[i]);


      output.writeToPixel(x, y, static_cast<char>((ray.t < 1e30f) * 255));
    }
  }

  output.save("./output.ppm");

  return 0;
}
