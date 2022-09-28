#pragma once
#include <stdint.h>

static uint32_t seed = 0x12345678;
uint32_t RandomUInt(){
  seed ^= seed << 13;
  seed ^= seed >> 17;
  seed ^= seed << 5;
  return seed;
}

float RandomFloat() { return RandomUInt() * 2.3283064365387e-10f; }
