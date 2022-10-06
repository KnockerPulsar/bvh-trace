#pragma once

#include "aabb.h"
#include "bvh.h"
#include "ray.h"
class BVHInstance {
	public:
		BVHInstance() = default;
		BVHInstance(bvt::BVH* blas) : bvh(blas) { SetTransform(mat4()); }
		void SetTransform(const mat4& transform);
		void Intersect(bvt::Ray& ray);

	private:
		bvt::BVH* bvh = nullptr;
		mat4 invTransform;
	public:
		AABB bounds;
};
