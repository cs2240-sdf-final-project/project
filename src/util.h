#ifndef UTIL_H
#define UTIL_H

#include "linmath.h"

extern int enzyme_dup;
extern int enzyme_dupnoneed;
extern int enzyme_out;
extern int enzyme_const;

static inline void ray_step(vec3 out, const vec3 origin, const vec3 direction, float t) {
    vec3 step_size;
    vec3_scale(step_size, direction, t);
    vec3_add(out, origin, step_size);
}

static inline float clamp(float x, float min, float max) {
    return fmaxf(fminf(x, max), min);
}

static inline void dehomogenize(vec3 out, const vec4 in) {
    for (int i = 0; i < 3; i++) {
        out[i] = in[i] / in[3];
    }
}

static inline float lerp(float x, float in_min, float in_max, float out_min, float out_max) {
    return out_min + (x - in_min) * (out_max - out_min) / (in_max - in_min);
}

static inline float vec3_distance(const vec3 a, const vec3 b) {
    vec3 displacement;
    vec3_sub(displacement, b, a);
    return vec3_len(displacement);
}

static inline void vec3_set(vec3 out, float value) {
    vec3 to_set = { value, value, value };
    vec3_dup(out, to_set);
}

static inline void vec3_componentwise_mul(vec3 out, const vec3 a, const vec3 b) {
    for (int i = 0; i < 3; i++) {
        out[i] = a[i] * b[i];
    }
}

static inline void vec3_set(vec3 out, float v1, float v2, float v3) {
    vec3 to_set = { v1, v2, v3 };
    vec3_dup(out, to_set);
}

static inline void vec2_abs(vec2 out, const vec2 in) {
    for (int i = 0; i < 2; i++) {
        out[i] = fabsf(in[i]);
    }
}

static inline void vec3_cwiseProduct(vec3 out, vec3 a, vec3 b){
    out[0] = a[0]*b[0];
    out[1] = a[1]*b[1];
    out[2] = a[2]*b[2];
}

static inline long lclamp(long x, long min, long max) {
    if (x > max) {
        return max;
    }
    if (x < min) {
        return min;
    }
    return x;
}

static inline float dot2(vec3 v) {
    return vec3_mul_inner(v, v);
}

static inline float sign(float x) {
    return (x > 0.0f) - (x < 0.0f);
}

#endif // UTIL_H
