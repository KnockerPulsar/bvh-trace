#pragma once
#include "vec_math.h"

struct AABB {
	float3 min = make_float3(1e30f), max = make_float3(-1e30f);
	void grow(const float3& p) { min = fminf(min, p), max = fmaxf(max, p); }	

	void grow(const AABB& other) { 
		if(other.min.x != 1e30f)
			grow(other.min), grow(other.max);
	}	

	float area() const {
		float3 e = max - min;
		return e.x * e.y + e.x * e.z + e.y * e.z;
	}
};
