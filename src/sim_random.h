#ifndef SIM_RANDOM_H
#define SIM_RANDOM_H

#include <stdint.h>
#include "linmath.h"

typedef struct {
    uint64_t s[1];
} RandomState;

// https://nullprogram.com/blog/2017/09/21/
uint32_t spcg32(uint64_t s[1]);

RandomState *make_random(void);
void free_random(RandomState *rng);

// Returns a random float between 0 and 1.
float random_next_float(RandomState *state);

long sample_binomial(float prob, RandomState *rng);

void sample_hemisphere(vec3 out, const vec3 normal, RandomState *rng);

float generate_normal_random(RandomState* random);

#endif
