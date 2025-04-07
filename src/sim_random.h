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
    uint64_t shift = 29 - (*s >> 61);
    uint64_t ret = *s >> shift;
    return (uint32_t)ret;
}

typedef struct {
    uint64_t s[1];
} RandomState;

RandomState make_random(void) {
    RandomState state;
    uint64_t s = UINT64_C(16018692385596851550);
    state.s[0] = s;
    return state;
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
