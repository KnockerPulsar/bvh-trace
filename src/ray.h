#pragma once
#include "vec_math.h"

struct Ray { float3 O, D, rD; float t = 1e30f; };
