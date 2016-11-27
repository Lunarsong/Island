#pragma once

#include <stdint.h>
#include <vector>

#include "noise/simplex_noise.h"

struct IslandParams {
  uint32_t seed = 0;
  uint32_t size_x;
  uint32_t size_y;

  float min_height;
  float max_height;

  int num_height_maps = 1;
};

struct IslandGenerator {
  IslandParams island_params;

  std::vector<SimplexNoise> simplex_gens;

  static IslandGenerator FromParams(const IslandParams& params);

  float HeightCentered(float x, float y) const;
  float Height(float x, float y) const;
  float Height(uint32_t x, uint32_t y) const;
};

float* GenerateHeightmapRegion(const IslandGenerator& gen, uint32_t start_x,
                               uint32_t start_y, uint32_t size_x,
                               uint32_t size_y, uint32_t intervals = 1);
float* GenerateHeightmap(uint32_t seed, uint32_t size_x, uint32_t size_y);
float* GenerateIslandHeightmap(uint32_t seed, uint32_t size_x, uint32_t size_y);
float* GenerateIslandHeightMask(uint32_t seed, uint32_t size_x,
                                uint32_t size_y);
