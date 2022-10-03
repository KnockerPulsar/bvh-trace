#pragma once

#include "defs.h"
#include "bvh_node.h"
#include "tri.h"

#define N 12583

void UpdateNodeBounds(uint nodeIndex);
void Subdivide(uint nodeIndex);
void BuildBVH();
