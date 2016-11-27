#pragma once

class SimplexNoise {
 public:
  SimplexNoise();
  ~SimplexNoise();

  void SetSeed(unsigned int seed);
  void SetBounds(float min, float max);

  void SetOctaves(float octaves);
  void SetPersistence(float persistence);
  void SetScale(float scale);

  unsigned int GetSeed() const;

  const void GetBounds(float& min, float& nax);

  float GetOctaves() const;
  float GetPersistence() const;
  float GetScale() const;

  // Noise functions
  float Noise(float x, float y) const;

 private:
  unsigned int seed_;

  // Return value range.
  float bounds_[2];

  // Noise settings.
  float octaves_;
  float persistence_;
  float scale_;

  int permutation_table_[512];

  float NoiseOctave2D(const float octaves, const float persistence,
                      const float scale, const float x, const float y) const;
  float NoiseRaw2D(const float x, const float y) const;
};
