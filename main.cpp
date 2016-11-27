#include "island/island.h"

#include <cmath>
#include <iostream>

int main() {
  const uint32_t size = pow(2, 13) + 1;

  IslandParams params;
  params.size_x = size;
  params.size_y = size;
  params.num_height_maps = 1;
  params.seed = time(NULL);
  params.min_height = -100.0f;
  params.max_height = 256.0f;
  IslandGenerator island_gen = IslandGenerator::FromParams(params);

  std::cout << "Island seed: " << params.seed << std::endl;
  std::cout << "Island size: " << size << std::endl;

  float* heightmap = GenerateHeightmapRegion(island_gen, 0, 0, size, size, 1);

  delete[] heightmap;

  return 0;
}