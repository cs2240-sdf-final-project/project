#include <stdint.h>
#include <stdio.h>
#include "linmath.h"

#include "sim_random.h"

// https://nullprogram.com/blog/2017/09/21/
uint32_t spcg32(uint64_t s[1])
{
    uint64_t m = 0x9b60933458e17d7d;
    uint64_t a = 0xd737232eeccdf7ed;
    *s = *s * m + a;
    int shift = 29 - (*s >> 61);
    uint32_t ret = (uint32_t)((*s) >> shift);
    return ret;
}

void make_random(RandomState *ret) {
    uint64_t s = UINT64_C(16018692385596851550);
    ret->s[0] = s;
}

// Returns a random float between 0 and 1.
float random_next_float(RandomState *state) {
    uint32_t scale = spcg32(state->s);
    double value = (double)scale;
    double max = (double)UINT32_MAX;
    double ret = value / max;
    return (float)ret;
}

long sample_binomial(float prob, RandomState *rng) {
    long ret = 0;

    for (;;) {
        float sample = random_next_float(rng);
        if (sample < prob) {
            ret += 1;
        } else {
            break;
        }
    }
    return ret;
}

////////// BEGIN CODE FROM CHATGPT
/* Marsaglia (1972) algorithm: uniform point on a unit sphere */
static void
sample_unit_sphere(vec3 out, RandomState *rng)
{
    for (;;) {
        float x1 = random_next_float(rng) * 2.0f - 1.0f;   /* (-1, 1) */
        float x2 = random_next_float(rng) * 2.0f - 1.0f;
        float sum = x1 * x1 + x2 * x2;
        if (sum >= 1.0f) continue;                         /* outside unit disc */

        float factor = 2.0f * sqrtf(1.0f - sum);
        out[0] = x1 * factor;
        out[1] = x2 * factor;
        out[2] = 1.0f - 2.0f * sum;
        return;                                           /* ‖out‖₂ == 1 */
    }
}
void
sample_hemisphere(vec3 out, const vec3 normal, RandomState *rng)
{
    /* 1. sample whole sphere */
    sample_unit_sphere(out, rng);

    if (vec3_mul_inner(out, normal) < 0.0f) {
        vec3_scale(out, out, -1.0f);       /* mirror through the origin */
    }
}
////////// END CODE FROM CHATGPT
