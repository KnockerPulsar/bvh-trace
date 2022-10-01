#include "vec_math.h"

struct AABB {
	float3 min = make_float3(1e30f), max = make_float3(-1e30f);
	void grow(float3 p) { min = fminf(min, p), max = fmaxf(max, p); }	
	float area() const {
		float3 e = max - min;
		return e.x * e.y + e.x * e.z + e.y * e.z;
	}
};
