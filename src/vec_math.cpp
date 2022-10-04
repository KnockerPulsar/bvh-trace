#include "vec_math.h"

float3 make_float3(float a) {
  return float3 { a, a, a };  
}

float3 make_float3(float a, float b, float c) {
  return float3 { a, b, c };  
}

float3 make_float3(float a, float b, float c, float d) {
  return float3 { a, b, c, d };  
}

float3 make_float3(float3 a, float w) {
  return float3 { a.x, a.y, a.z, w };  
}

float3 operator*(const float3& a, float b) {
  return float3{a.x * b, a.y * b, a.z * b};
}

float3 operator*(const float3& a, const float3& b) {
  return float3{a.x * b.x, a.y * b.y, a.z * b.z};
}

float3 operator+(const float3& a, const float3& b) {
  return float3{a.x + b.x, a.y + b.y, a.z + b.z};
}

float3 operator-(const float3& a, const float3& b) {
  return float3{a.x - b.x, a.y - b.y, a.z - b.z};
}

float3 cross(const float3& a, const float3& b) {
  return float3{
    a.y * b.z - a.z * b.y,
      -(a.x * b.z - b.x * a.z),
      a.x * b.y - b.x * a.y
  };
}

float dot(const float3& a, const float3& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

float3 normalize(const float3 a) {
  float3 invLen = make_float3(1 / sqrtf(dot(a,a)));
  return a * invLen;
}

float3 fminf(const float3& a, const float3& b) {
  return float3 {
    fminf(a.x, b.x),
      fminf(a.y, b.y),
      fminf(a.z, b.z)
  };
}

float3 fmaxf(const float3& a, const float3& b) {
  return float3 { 
    fmaxf(a.x, b.x),
      fmaxf(a.y, b.y),
      fmaxf(a.z, b.z)
  };
}


mat4 mat4::Translate(const float x, const float y, const float z) {
  mat4 r;
  r.cell[3] = x;
  r.cell[7] = y;
  r.cell[11] = z;

  return r;
}

mat4 mat4::Translate( const float3& p) {
  return mat4::Translate(p.x, p.y, p.z);
}

mat4 mat4::RotateY(float a) {
  return {
    cosf(a), 0, -sinf(a), 0,
      0      , 1, 0       , 0,
      sinf(a), 0, cosf(a) , 0,
      0      , 0, 0       , 1
  };
}

mat4 mat4::Inverted() const {
  double inv[16], det;
  int i;

  inv[0] = cell[5]  * cell[10] * cell[15] - 
    cell[5]  * cell[11] * cell[14] - 
    cell[9]  * cell[6]  * cell[15] + 
    cell[9]  * cell[7]  * cell[14] +
    cell[13] * cell[6]  * cell[11] - 
    cell[13] * cell[7]  * cell[10];

  inv[4] = -cell[4]  * cell[10] * cell[15] + 
    cell[4]  * cell[11] * cell[14] + 
    cell[8]  * cell[6]  * cell[15] - 
    cell[8]  * cell[7]  * cell[14] - 
    cell[12] * cell[6]  * cell[11] + 
    cell[12] * cell[7]  * cell[10];

  inv[8] = cell[4]  * cell[9] * cell[15] - 
    cell[4]  * cell[11] * cell[13] - 
    cell[8]  * cell[5] * cell[15] + 
    cell[8]  * cell[7] * cell[13] + 
    cell[12] * cell[5] * cell[11] - 
    cell[12] * cell[7] * cell[9];

  inv[12] = -cell[4]  * cell[9] * cell[14] + 
    cell[4]  * cell[10] * cell[13] +
    cell[8]  * cell[5] * cell[14] - 
    cell[8]  * cell[6] * cell[13] - 
    cell[12] * cell[5] * cell[10] + 
    cell[12] * cell[6] * cell[9];

  inv[1] = -cell[1]  * cell[10] * cell[15] + 
    cell[1]  * cell[11] * cell[14] + 
    cell[9]  * cell[2] * cell[15] - 
    cell[9]  * cell[3] * cell[14] - 
    cell[13] * cell[2] * cell[11] + 
    cell[13] * cell[3] * cell[10];

  inv[5] = cell[0]  * cell[10] * cell[15] - 
    cell[0]  * cell[11] * cell[14] - 
    cell[8]  * cell[2] * cell[15] + 
    cell[8]  * cell[3] * cell[14] + 
    cell[12] * cell[2] * cell[11] - 
    cell[12] * cell[3] * cell[10];

  inv[9] = -cell[0]  * cell[9] * cell[15] + 
    cell[0]  * cell[11] * cell[13] + 
    cell[8]  * cell[1] * cell[15] - 
    cell[8]  * cell[3] * cell[13] - 
    cell[12] * cell[1] * cell[11] + 
    cell[12] * cell[3] * cell[9];

  inv[13] = cell[0]  * cell[9] * cell[14] - 
    cell[0]  * cell[10] * cell[13] - 
    cell[8]  * cell[1] * cell[14] + 
    cell[8]  * cell[2] * cell[13] + 
    cell[12] * cell[1] * cell[10] - 
    cell[12] * cell[2] * cell[9];

  inv[2] = cell[1]  * cell[6] * cell[15] - 
    cell[1]  * cell[7] * cell[14] - 
    cell[5]  * cell[2] * cell[15] + 
    cell[5]  * cell[3] * cell[14] + 
    cell[13] * cell[2] * cell[7] - 
    cell[13] * cell[3] * cell[6];

  inv[6] = -cell[0]  * cell[6] * cell[15] + 
    cell[0]  * cell[7] * cell[14] + 
    cell[4]  * cell[2] * cell[15] - 
    cell[4]  * cell[3] * cell[14] - 
    cell[12] * cell[2] * cell[7] + 
    cell[12] * cell[3] * cell[6];

  inv[10] = cell[0]  * cell[5] * cell[15] - 
    cell[0]  * cell[7] * cell[13] - 
    cell[4]  * cell[1] * cell[15] + 
    cell[4]  * cell[3] * cell[13] + 
    cell[12] * cell[1] * cell[7] - 
    cell[12] * cell[3] * cell[5];

  inv[14] = -cell[0]  * cell[5] * cell[14] + 
    cell[0]  * cell[6] * cell[13] + 
    cell[4]  * cell[1] * cell[14] - 
    cell[4]  * cell[2] * cell[13] - 
    cell[12] * cell[1] * cell[6] + 
    cell[12] * cell[2] * cell[5];

  inv[3] = -cell[1] * cell[6] * cell[11] + 
    cell[1] * cell[7] * cell[10] + 
    cell[5] * cell[2] * cell[11] - 
    cell[5] * cell[3] * cell[10] - 
    cell[9] * cell[2] * cell[7] + 
    cell[9] * cell[3] * cell[6];

  inv[7] = cell[0] * cell[6] * cell[11] - 
    cell[0] * cell[7] * cell[10] - 
    cell[4] * cell[2] * cell[11] + 
    cell[4] * cell[3] * cell[10] + 
    cell[8] * cell[2] * cell[7] - 
    cell[8] * cell[3] * cell[6];

  inv[11] = -cell[0] * cell[5] * cell[11] + 
    cell[0] * cell[7] * cell[9] + 
    cell[4] * cell[1] * cell[11] - 
    cell[4] * cell[3] * cell[9] - 
    cell[8] * cell[1] * cell[7] + 
    cell[8] * cell[3] * cell[5];

  inv[15] = cell[0] * cell[5] * cell[10] - 
    cell[0] * cell[6] * cell[9] - 
    cell[4] * cell[1] * cell[10] + 
    cell[4] * cell[2] * cell[9] + 
    cell[8] * cell[1] * cell[6] - 
    cell[8] * cell[2] * cell[5];

  det = cell[0] * inv[0] + cell[1] * inv[4] + cell[2] * inv[8] + cell[3] * inv[12];

  mat4 retVal;

  if(det != 0) {
    const float invDet = 1.0f / det;
    for (int i = 0; i < 16 ; i++) {
      retVal.cell[i] = inv[i] * invDet;
    }
  }

  return retVal;
}


float3 operator* (const float3& b, const mat4& a) {
  return make_float3( 
      a.cell[0] * b.x + a.cell[1] * b.y + a.cell[2] * b.z + a.cell[3] * b.w,
      a.cell[4] * b.x + a.cell[5] * b.y + a.cell[6] * b.z + a.cell[7] * b.w,
      a.cell[8] * b.x + a.cell[9] * b.y + a.cell[10] * b.z + a.cell[11] * b.w,
      a.cell[12] * b.x + a.cell[13] * b.y + a.cell[14] * b.z + a.cell[15] * b.w 
      );
}

mat4 operator* (const mat4& a, const mat4& b) {
  mat4 r;
  for (uint i =0; i < 16; i+=4) {
    for (uint j = 0; j < 4; ++j) {
      r[i+j] = 
        (a.cell[i + 0] * b.cell[j + 0]) + 
        (a.cell[i + 1] * b.cell[j + 4]) + 
        (a.cell[i + 2] * b.cell[j + 8]) + 
        (a.cell[i + 3] * b.cell[j + 12]);
    } 

  }

  return r;
}

float3 TransformPosition( const float3& a, const mat4& M) {
  return make_float3(a, 1) * M;
}

float3 TransformVector( const float3& a, const mat4& M) {
  return make_float3(a, 0) * M;
}
