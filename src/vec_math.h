#pragma once
#include <algorithm>
#include <numeric>
#include <math.h>
#define ALIGN(x) __attribute__( ( aligned(x) ))

struct mat4;

struct ALIGN(16) float3 {
  float3() = default;
  
  union { struct { float x,y,z, w; }; float cell[4]; };
  float operator[] (const int n) const { return cell[n]; }
};

float3 make_float3(float a); 
float3 make_float3(float a, float b, float c);
float3 make_float3(float a, float b, float c, float d);
float3 make_float3(float3 a, float w);

float3 operator*(const float3& a, float b);
float3 operator-(const float3& a, float b); 

float3 operator*(const float3& a, const float3& b);
float3 operator+(const float3& a, const float3& b);
float3 operator-(const float3& a, const float3& b);

float3 operator*=(float3& a, float b);
float3 operator+=(float3& a, const float3& b);

float3 cross(const float3& a, const float3& b);
float dot(const float3& a, const float3& b);
float3 normalize(const float3 a);

float3 fminf(const float3& a, const float3& b);
float3 fmaxf(const float3& a, const float3& b);

// Column major
struct ALIGN(16) mat4 {
  float cell[16] = {
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  };

  float& operator[] (const int idx) { return cell[idx]; }
  float operator() (const int i, const int j) const { return cell[i * 4 + j]; }
  float& operator() (const int i, const int j) { return cell[i * 4 + j]; }

  static mat4 Translate(const float x, const float y, const float z);
  static mat4 Translate( const float3& p);
  static mat4 RotateY(float a);
  static mat4 RotateX(float a);
  static mat4 RotateZ(float a);
  static mat4 Scale(float a);

  mat4 Inverted() const; 
};

float3 operator* (const float3& b, const mat4& a);
mat4 operator* (const mat4& a, const mat4& b);
float3 TransformPosition( const float3& a, const mat4& M);
float3 TransformVector( const float3& a, const mat4& M);
