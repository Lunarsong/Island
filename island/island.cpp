#include "island.h"

#include <cmath>
#include <cassert>

IslandGenerator IslandGenerator::FromParams(const IslandParams& params) {
  IslandGenerator generator;
  generator.island_params = params;

  for (int i = 0; i < params.num_height_maps; ++i) {
    SimplexNoise noise;
    noise.SetSeed(params.seed + i);
    noise.SetBounds(0.0f, 1.0f);
    noise.SetOctaves(4);
    noise.SetPersistence(0.8f);
    noise.SetScale(0.00025f);

    generator.simplex_gens.emplace_back(std::move(noise));
  }

  return generator;
}

float IslandGenerator::Height(uint32_t x, uint32_t y) const {
  // Radial mask params.
  const float half_size_x = (island_params.size_x - 1) / 2.0f;
  const float half_size_y = (island_params.size_y - 1) / 2.0f;

  const float min_size = std::min(half_size_x, half_size_y);
  const float island_radius = std::min(min_size * 0.9f, min_size - 2.0f);

  // Calculate radial mask value.
  const float dx = x - half_size_x;
  const float dy = y - half_size_y;
  const float distance = sqrt(dx * dx + dy * dy) / island_radius;

  const float island_mask = std::max(0.0f, 1.0f - distance);

  // Height map.
  float height_noise = 1.0f;
  for (auto& noise : simplex_gens) {
    height_noise *= noise.Noise(x, y);
  }

  const float height_range =
      island_params.max_height - island_params.min_height;
  const float final_height =
      height_noise * island_mask * height_range + island_params.min_height;

  return std::max(0.0f, final_height);
}

float IslandGenerator::HeightCentered(float x, float y) const {
  const float half_size_x = (island_params.size_x - 1) / 2.0f;
  const float half_size_y = (island_params.size_y - 1) / 2.0f;

  return Height(x + half_size_x, y + half_size_y);
}

float IslandGenerator::Height(float x, float y) const {
  // Radial mask params.
  const float half_size_x = (island_params.size_x - 1) / 2.0f;
  const float half_size_y = (island_params.size_y - 1) / 2.0f;

  const float min_size = std::min(half_size_x, half_size_y);
  const float island_radius = std::min(min_size * 0.9f, min_size - 2.0f);

  // Calculate radial mask value.
  const float dx = x - half_size_x;
  const float dy = y - half_size_y;
  const float distance = sqrt(dx * dx + dy * dy) / island_radius;

  const float island_mask = std::max(0.0f, 1.0f - distance);

  // Height map.
  float height_noise = 1.0f;
  for (auto& noise : simplex_gens) {
    height_noise *= noise.Noise(x, y);
  }

  const float height_range =
      island_params.max_height - island_params.min_height;
  const float final_height =
      height_noise * island_mask * height_range + island_params.min_height;

  return std::max(0.0f, final_height);
}

float* GenerateHeightmapRegion(const IslandGenerator& gen, uint32_t start_x,
                               uint32_t start_y, uint32_t size_x,
                               uint32_t size_y, uint32_t intervals) {
  assert(intervals > 0);

  float* heights = new float[size_x * size_y];

  for (uint32_t y = 0; y < size_y; ++y) {
    for (uint32_t x = 0; x < size_x; ++x) {
      const uint32_t array_index = y * size_x + x;
      heights[array_index] =
          gen.Height(start_x + x * intervals, start_y + y * intervals);
    }
  }

  return heights;
}

float* GenerateHeightmap(uint32_t seed, int num_octaves, uint32_t size_x,
                         uint32_t size_y) {
  assert(size_x);
  assert(size_y);

  // ExponentialDistributedPerlin exponential_distributed_perlin(seed);
  SimplexNoise noise;
  noise.SetSeed(seed);
  noise.SetBounds(0.0f, 1.0f);
  noise.SetOctaves(num_octaves);
  noise.SetScale(0.00025f);

  float* heights = new float[size_x * size_y];

  for (uint32_t y = 0; y < size_y; ++y) {
    for (uint32_t x = 0; x < size_x; ++x) {
      const uint32_t array_index = y * size_x + x;

      // Retrieve noise value from -1.0f to 1.0f.
      /*float noise = exponential_distributed_perlin.Noise2D(
          (float)x / (float)(size_x - 1), (float)y / (float)(size_y - 1),
          num_octaves);

      // Normalize the noise value.
      noise = noise * 0.5f + 0.5f;

      // Set the height.
      heights[array_index] = noise * height_range + min_height;
      // std::cout << noise << std::endl;*/

      heights[array_index] = noise.Noise(x, y);
    }
  }

  return heights;
}

float* GenerateIslandHeightmap(uint32_t seed, uint32_t size_x,
                               uint32_t size_y) {
  const float min_height = -20.0f;
  const float max_height = 256.0f;
  const float height_range = max_height - min_height;

  float* height_mask = GenerateIslandHeightMask(seed, size_x, size_y);
  float* heightmap = GenerateHeightmap(seed, 4, size_x, size_y);
  float* heightmap2 = GenerateHeightmap(seed + 1, 5, size_x, size_y);

  for (uint32_t y = 0; y < size_y; ++y) {
    for (uint32_t x = 0; x < size_x; ++x) {
      const uint32_t array_index = y * size_x + x;

      const float height_percent = heightmap[array_index] *
                                   heightmap2[array_index] *
                                   height_mask[array_index];
      const float height = height_percent * height_range + min_height;

      heightmap[array_index] = std::max(0.0f, height);
    }
  }

  delete[] height_mask;
  delete[] heightmap2;

  return heightmap;
}

float* GenerateIslandHeightMask(uint32_t seed, uint32_t size_x,
                                uint32_t size_y) {
  float* height_mask = new float[size_x * size_y];

  const float half_size_x = (size_x - 1) / 2.0f;
  const float half_size_y = (size_y - 1) / 2.0f;

  const float min_size = std::min(half_size_x, half_size_y);
  const float island_radius = std::min(min_size * 0.9f, min_size - 2.0f);

  for (uint32_t y = 0; y < size_y; ++y) {
    for (uint32_t x = 0; x < size_x; ++x) {
      const uint32_t array_index = y * size_x + x;

      float dx = x - half_size_x;
      float dy = y - half_size_y;

      float distance = sqrt(dx * dx + dy * dy) / island_radius;

      height_mask[array_index] = std::max(0.0f, 1.0f - distance);
    }
  }

  return height_mask;
}