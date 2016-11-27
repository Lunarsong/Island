#include <math.h>
#include <stdlib.h>
#include <ctime>
#include "simplex_noise.h"

// The gradients are the midpoints of the vertices of a cube.
static const int grad3[12][3] = {
    {1, 1, 0},  {-1, 1, 0},  {1, -1, 0}, {-1, -1, 0}, {1, 0, 1},  {-1, 0, 1},
    {1, 0, -1}, {-1, 0, -1}, {0, 1, 1},  {0, -1, 1},  {0, 1, -1}, {0, -1, -1}};

// The gradients are the midpoints of the vertices of a hypercube.
static const int grad4[32][4] = {
    {0, 1, 1, 1},  {0, 1, 1, -1},  {0, 1, -1, 1},  {0, 1, -1, -1},
    {0, -1, 1, 1}, {0, -1, 1, -1}, {0, -1, -1, 1}, {0, -1, -1, -1},
    {1, 0, 1, 1},  {1, 0, 1, -1},  {1, 0, -1, 1},  {1, 0, -1, -1},
    {-1, 0, 1, 1}, {-1, 0, 1, -1}, {-1, 0, -1, 1}, {-1, 0, -1, -1},
    {1, 1, 0, 1},  {1, 1, 0, -1},  {1, -1, 0, 1},  {1, -1, 0, -1},
    {-1, 1, 0, 1}, {-1, 1, 0, -1}, {-1, -1, 0, 1}, {-1, -1, 0, -1},
    {1, 1, 1, 0},  {1, 1, -1, 0},  {1, -1, 1, 0},  {1, -1, -1, 0},
    {-1, 1, 1, 0}, {-1, 1, -1, 0}, {-1, -1, 1, 0}, {-1, -1, -1, 0}};

// A lookup table to traverse the simplex around a given point in 4D.
static const int simplex[64][4] = {
    {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 0, 0, 0}, {0, 2, 3, 1}, {0, 0, 0, 0},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 2, 3, 0}, {0, 2, 1, 3}, {0, 0, 0, 0},
    {0, 3, 1, 2}, {0, 3, 2, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {1, 3, 2, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 2, 0, 3},
    {0, 0, 0, 0}, {1, 3, 0, 2}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {2, 3, 0, 1}, {2, 3, 1, 0}, {1, 0, 2, 3}, {1, 0, 3, 2}, {0, 0, 0, 0},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {2, 0, 3, 1}, {0, 0, 0, 0}, {2, 1, 3, 0},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {2, 0, 1, 3}, {0, 0, 0, 0},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {3, 0, 1, 2}, {3, 0, 2, 1}, {0, 0, 0, 0},
    {3, 1, 2, 0}, {2, 1, 0, 3}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {3, 1, 0, 2}, {0, 0, 0, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}};

int fastfloor(const float x) { return x > 0 ? (int)x : (int)x - 1; }

float dot(const int* g, const float x, const float y) {
  return g[0] * x + g[1] * y;
}

float dot(const int* g, const float x, const float y, const float z) {
  return g[0] * x + g[1] * y + g[2] * z;
}

float dot(const int* g, const float x, const float y, const float z,
          const float w) {
  return g[0] * x + g[1] * y + g[2] * z + g[3] * w;
}

float random(float max) {
  int r;
  float s;

  r = rand();
  s = (float)(r & 0x7fff) / (float)0x7fff;

  return (s * max);
}

SimplexNoise::SimplexNoise() {
  SetSeed(time(NULL));

  bounds_[0] = -1.0f;
  bounds_[1] = 1.0f;

  scale_ = 1.0f;
  persistence_ = 1.0f;
  octaves_ = 0.0f;
}

SimplexNoise::~SimplexNoise() {}

void SimplexNoise::SetSeed(unsigned int seed) {
  seed_ = seed;
  srand(seed_);

  // Set up the random numbers table
  for (int i = 0; i < 256; ++i) {
    permutation_table_[i] = i;  // put each number in once
  }

  // Randomize the random numbers table
  for (int i = 0; i < 256; ++i) {
    int j = (int)random(256);
    int k = permutation_table_[i];
    permutation_table_[i] = permutation_table_[j];
    permutation_table_[j] = k;

    permutation_table_[256 + i] = permutation_table_[j];
    permutation_table_[256 + j] = k;
  }
}

void SimplexNoise::SetBounds(float min, float max) {
  bounds_[0] = min;
  bounds_[1] = max;
}

void SimplexNoise::SetOctaves(float octaves) { octaves_ = octaves; }

void SimplexNoise::SetPersistence(float persistence) {
  persistence_ = persistence;
}

void SimplexNoise::SetScale(float scale) { scale_ = scale; }

unsigned int SimplexNoise::GetSeed() const { return seed_; }

const void SimplexNoise::GetBounds(float& min, float& max) {
  min = bounds_[0];
  max = bounds_[1];
}

float SimplexNoise::GetOctaves() const { return octaves_; }

float SimplexNoise::GetPersistence() const { return persistence_; }

float SimplexNoise::GetScale() const { return scale_; }

// Noise functions
float SimplexNoise::Noise(float x, float y) const {
  return NoiseOctave2D(octaves_, persistence_, scale_, x, y) *
             (bounds_[1] - bounds_[0]) * .5f +
         (bounds_[1] + bounds_[0]) * 0.5f;
}

float SimplexNoise::NoiseRaw2D(const float x, const float y) const {
  // Noise contributions from the three corners
  float n0, n1, n2;

  // Skew the input space to determine which simplex cell we're in
  float F2 = 0.5 * (sqrtf(3.0) - 1.0);
  // Hairy factor for 2D
  float s = (x + y) * F2;
  int i = fastfloor(x + s);
  int j = fastfloor(y + s);

  float G2 = (3.0 - sqrtf(3.0)) / 6.0;
  float t = (i + j) * G2;
  // Unskew the cell origin back to (x,y) space
  float X0 = i - t;
  float Y0 = j - t;
  // The x,y distances from the cell origin
  float x0 = x - X0;
  float y0 = y - Y0;

  // For the 2D case, the simplex shape is an equilateral triangle.
  // Determine which simplex we are in.
  int i1, j1;  // Offsets for second (middle) corner of simplex in (i,j) coords
  if (x0 > y0) {
    i1 = 1;
    j1 = 0;
  }  // lower triangle, XY order: (0,0)->(1,0)->(1,1)
  else {
    i1 = 0;
    j1 = 1;
  }  // upper triangle, YX order: (0,0)->(0,1)->(1,1)

  // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
  // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
  // c = (3-sqrt(3))/6
  float x1 =
      x0 - i1 + G2;  // Offsets for middle corner in (x,y) unskewed coords
  float y1 = y0 - j1 + G2;
  float x2 =
      x0 - 1.0 + 2.0 * G2;  // Offsets for last corner in (x,y) unskewed coords
  float y2 = y0 - 1.0 + 2.0 * G2;

  // Work out the hashed gradient indices of the three simplex corners
  int ii = i & 255;
  int jj = j & 255;
  int gi0 = permutation_table_[ii + permutation_table_[jj]] % 12;
  int gi1 = permutation_table_[ii + i1 + permutation_table_[jj + j1]] % 12;
  int gi2 = permutation_table_[ii + 1 + permutation_table_[jj + 1]] % 12;

  // Calculate the contribution from the three corners
  float t0 = 0.5 - x0 * x0 - y0 * y0;
  if (t0 < 0)
    n0 = 0.0;
  else {
    t0 *= t0;
    n0 = t0 * t0 *
         dot(grad3[gi0], x0, y0);  // (x,y) of grad3 used for 2D gradient
  }

  float t1 = 0.5 - x1 * x1 - y1 * y1;
  if (t1 < 0)
    n1 = 0.0;
  else {
    t1 *= t1;
    n1 = t1 * t1 * dot(grad3[gi1], x1, y1);
  }

  float t2 = 0.5 - x2 * x2 - y2 * y2;
  if (t2 < 0)
    n2 = 0.0;
  else {
    t2 *= t2;
    n2 = t2 * t2 * dot(grad3[gi2], x2, y2);
  }

  // Add contributions from each corner to get the final noise value.
  // The result is scaled to return values in the interval [-1,1].
  return 70.0 * (n0 + n1 + n2);
}

float SimplexNoise::NoiseOctave2D(const float octaves, const float persistence,
                                  const float scale, const float x,
                                  const float y) const {
  float total = 0;
  float frequency = scale;
  float amplitude = 1;

  // We have to keep track of the largest possible amplitude,
  // because each octave adds more, and we need a value in [-1, 1].
  float maxAmplitude = 0;

  for (int i = 0; i < octaves; i++) {
    total += NoiseRaw2D(x * frequency, y * frequency) * amplitude;

    frequency *= 2;
    maxAmplitude += amplitude;
    amplitude *= persistence;
  }

  return total / maxAmplitude;
}
