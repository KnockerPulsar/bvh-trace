// https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/
#include "bvh_instance.h"
#include "vec_math.h"
#include "random.h"
#include "PPMImage.h"
#include "bvh_node.h"
#include "bvh.h"
#include "tri.h"
#include "ray.h"
#include "tlas.h"
#include "config.h"
#include "util.h"

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

#define INSTANCE_COUNT 128
extern BVHInstance bvhInstance[INSTANCE_COUNT];
extern TLAS tlas;

float3 *position, *direction, *orientation;

void Animate() {
  for (int i = 0 ; i < INSTANCE_COUNT; i++) {
    mat4 R = mat4::RotateX(orientation[i].x) 
      * mat4::RotateY(orientation[i].y)
      * mat4::RotateZ(orientation[i].z)
      * mat4::Scale(0.2f);

    bvhInstance[i].SetTransform(mat4::Translate(position[i]) * R);
    position[i] += direction[i], orientation[i] += direction[i];;

    if(position[i].x < -3 || position[i].x > 3) direction[i].x *= -1;
    if(position[i].y < -3 || position[i].y > 3) direction[i].y *= -1;
    if(position[i].z < -3 || position[i].z > 3) direction[i].z *= -1;
  }
}

int main() {

  bvt::Image output = bvt::Image::fromColor(RES_X, RES_Y, make_float3(0));
  bvt::BVH* bvh = new bvt::BVH("assets/armadillo.tri", 30000);

  for (int i = 0; i < INSTANCE_COUNT; i++) {
    bvhInstance[i] = BVHInstance(bvh);
  }

  tlas = TLAS(bvhInstance, INSTANCE_COUNT);
  position = new float3[INSTANCE_COUNT];
  direction = new float3[INSTANCE_COUNT];
  orientation = new float3[INSTANCE_COUNT];

  for (int i = 0 ; i < INSTANCE_COUNT; i++) {
    position[i] = make_float3(RandomFloat(), RandomFloat(), RandomFloat()) - 0.5f;
    position[i] *= 4;

    direction[i] = normalize(position[i]) * 0.05f;
    orientation[i] = make_float3(RandomFloat(), RandomFloat(), RandomFloat()) * 2.5f;
  }


  InitWindow(RES_X, RES_Y, "BVH Tracer");

  Texture frameBuffer = LoadTextureFromImage({
      .data = (void*)output.getData().data(),             
      .width = (int) RES_X,              
      .height = (int) RES_Y,
      .mipmaps = 1,
      .format = RL_PIXELFORMAT_UNCOMPRESSED_R8G8B8
      });


  const int tile_size = 8;

  float3 posOffset = {0};
  float angleX = 0.0f, angleY = 0;

  while(!WindowShouldClose()) {

    Animate();
    tlas.Build();

    if(IsKeyDown(KEY_D)) posOffset.x += 0.1f;
    if(IsKeyDown(KEY_A)) posOffset.x -= 0.1f;
    
    if(IsKeyDown(KEY_W)) posOffset.z += 0.1f;
    if(IsKeyDown(KEY_S)) posOffset.z -= 0.1f;
    
    if(IsKeyDown(KEY_E)) posOffset.y += 0.1f;
    if(IsKeyDown(KEY_Q)) posOffset.y -= 0.1f;

    if(IsKeyDown(KEY_UP)) angleX += 0.1f;
    if(IsKeyDown(KEY_DOWN)) angleX -= 0.1f;

    if(IsKeyDown(KEY_LEFT)) angleY += 0.1f;
    if(IsKeyDown(KEY_RIGHT)) angleY -= 0.1f;


    float3 p0 = TransformPosition(make_float3(-1,+1,2), mat4::RotateX(angleX) * mat4::RotateY(angleY));
    float3 p1 = TransformPosition(make_float3(+1,+1,2), mat4::RotateX(angleX) * mat4::RotateY(angleY));
    float3 p2 = TransformPosition(make_float3(-1,-1,2), mat4::RotateX(angleX) * mat4::RotateY(angleY));

    /* printf("Ray origin: (%f, %f, %f)\n", ray.O.x, ray.O.y, ray.O.z); */
    /* printf("Angle X: %f, Angle Y: %f\n", angleX, angleY); */


    auto start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for schedule(dynamic)
    for (int tile = 0; tile < (RES_X * RES_Y) / 64; tile++) {

      int x = tile % (RES_X / 8), y = tile / (RES_X / 8);

      bvt::Ray ray;
      ray.O = make_float3(0,0, -6.5) + posOffset;

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

        float3 pixelPos = ray.O + p0 
          + (p1-p0) * ((x * 8 + u)/(float)RES_X) 
          + (p2-p0) * ((y * 8 + v)/(float)RES_Y);

        ray.D = normalize(pixelPos - ray.O);
        ray.t = 1e30f;

        ray.rD = make_float3(1/ray.D.x , 1/ray.D.y, 1/ray.D.z); 

        tlas.Intersect(ray);

        uint c = ray.t < 1e30f ? (int)(255 / (1 + std::max( 0.f, ray.t - 4 ))) : 0;

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

