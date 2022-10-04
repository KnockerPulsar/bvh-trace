// https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/
#include "vec_math.h"
#include "random.h"
#include "PPMImage.h"
#include "bvh_node.h"
#include "bvh.h"
#include "tri.h"
#include "ray.h"
#include "tlas.h"
#include "config.h"

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

extern bvt::BVH bvh[2];
extern TLAS tlas;


/* void Animate() { */
/*   static float r = 0; */
/*   if((r+=0.05f) > 2 * PI) r -= 2 * PI; */
/*  */
/*   float a = sinf(r) * 0.5f; */
/*   for (int i = 0; i < N; i++) for (int j = 0; j < 3; j++) { */
/*  */
/*     // So you can access vertex0, vertex1, and vertex2 */
/*     // with iteration instead of having to write them manually */
/*     float3 o = (&original[i].vertex0)[j]; */
/*  */
/*     // Makes the rotation angle scale with the y axis position (height) */
/*     // of the vertex so the base doesn't move but the edges move significantly. */
/*     float s = a * (o.y - 0.2f) * 0.2f; */
/*  */
/*     // Rotating about the z axis? */
/*     float x = o.x * cosf( s ) - o.y * sinf( s ); */
/*     float y = o.x * sinf( s ) + o.y * cosf( s ); */
/*  */
/*     (&tri[i].vertex0)[j] = make_float3(x,y,o.z); */
/*   } */
/* } */

int main() {

  bvt::Image output = bvt::Image::fromColor(RES_X, RES_Y, make_float3(0));

  float3 p0{-1, 1, 2}, p1{1, 1, 2}, p2{-1, -1, 2};
  bvh[0] = bvt::BVH("assets/bigben.tri", 20944);
  bvh[1] = bvt::BVH("assets/bigben.tri", 20944);

  tlas = TLAS(bvh, 2);
  tlas.Build();

  InitWindow(RES_X, RES_Y, "BVH Tracer");

  Texture frameBuffer = LoadTextureFromImage({
      .data = (void*)output.getData().data(),             
      .width = (int) RES_X,              
      .height = (int) RES_Y,
      .mipmaps = 1,
      .format = RL_PIXELFORMAT_UNCOMPRESSED_R8G8B8
      });


  const int tile_size = 8;

  while(!WindowShouldClose()) {
    static float angle = 0;
    angle += 0.01f; if (angle > 2 * PI) angle -= 2 * PI;

    bvh[0].SetTranform(mat4::Translate(-1.3f, 0, 0));
    bvh[1].SetTranform(mat4::Translate(1.3f, 0, 0) * mat4::RotateY(angle));

#pragma parallel for schedule(dynamic)

    auto start = std::chrono::high_resolution_clock::now();
    for (int tile = 0; tile < 6400; tile++) {

      int x = tile % 80, y = tile / 80;
      bvt::Ray ray;
      ray.O = make_float3(1, 3.5f, -4.5f);

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

        tlas.Intersect(ray);

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

