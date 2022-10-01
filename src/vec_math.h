#pragma once
#include <numeric>
#include <math.h>
#define ALIGN(x) __attribute__( ( aligned(x) ))

struct ALIGN(16) float3 {
  float3() = default;
  
  union { struct { float x,y,z, dummy; }; float cell[4]; };
  float operator[] (const int n) const { return cell[n]; }
};

inline float3 make_float3(float a) {
  return float3 { a, a, a };  
}

inline float3 make_float3(float a, float b, float c) {
  return float3 { a, b, c };  
}

inline float3 operator*(const float3& a, float b) {
  return float3{a.x * b, a.y * b, a.z * b};
}

inline float3 operator*(const float3& a, const float3& b) {
  return float3{a.x * b.x, a.y * b.y, a.z * b.z};
}

inline float3 operator+(const float3& a, const float3& b) {
  return float3{a.x + b.x, a.y + b.y, a.z + b.z};
}

inline float3 operator-(const float3& a, const float3& b) {
  return float3{a.x - b.x, a.y - b.y, a.z - b.z};
}

inline float3 cross(const float3& a, const float3& b) {
  return float3{
      a.y * b.z - a.z * b.y,
      -(a.x * b.z - b.x * a.z),
      a.x * b.y - b.x * a.y
  };
}

inline float dot(const float3& a, const float3& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline float3 normalize(const float3 a) {
  float3 invLen = make_float3(1 / sqrtf(dot(a,a)));
  return a * invLen;
}

inline float3 fminf(const float3& a, const float3& b) {
  return float3 {
    fminf(a.x, b.x),
    fminf(a.y, b.y),
    fminf(a.z, b.z)
  };
}

inline float3 fmaxf(const float3& a, const float3& b) {
  return float3 { 
    fmaxf(a.x, b.x),
    fmaxf(a.y, b.y),
    fmaxf(a.z, b.z)
  };
}
