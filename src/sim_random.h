#ifndef SIM_RANDOM_H
#define SIM_RANDOM_H

#include <stdint.h>
#include <limits.h>

// https://nullprogram.com/blog/2017/09/21/
uint32_t spcg32(uint64_t s[1])
{
    uint64_t m = 0x9b60933458e17d7d;
    uint64_t a = 0xd737232eeccdf7ed;
    *s = *s * m + a;
    int shift = 29 - (*s >> 61);
    return *s >> shift;
}

typedef struct {
    uint64_t s[1];
} RandomState;

void random_init(RandomState *state) {
    uint64_t s = UINT64_C(16018692385596851550);
    state->s[0] = s;
}

// Returns a random float between 0 and 1.
float random_next_float(RandomState *state) {
    uint32_t scale = spcg32(state->s);
    double value = (double)scale;
    double max = (double)UINT32_MAX;
    double ret = value / max;
    return (float)ret;
}


#endif
