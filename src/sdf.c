#include <math.h>
#include "linmath.h"

float sdfCylinder(const vec3 pos,float radius, float height) {
    vec2 xz;
    xz[0] = pos[0];
    xz[1] = pos[2];
    float xzLen = vec2_len(xz);

    vec2 d;
    d[0] = xzLen - radius;
    d[1] = fabsf(pos[1]) - height;

    vec2 abs_d;
    vec2_abs(abs_d, d);

    float max_d = fmaxf(d[0], d[1]);
    float min_d = fminf(max_d, 0.0f);

    vec2 d_clamped;
    d_clamped[0] = fmaxf(d[0], 0.0f);
    d_clamped[1] = fmaxf(d[1], 0.0f);
    float dist = min_d + vec2_len(d_clamped);

    return dist;
}

float sdfPlane(const vec3 pos, const vec3 normal, float height) {
    float dist = vec3_mul_inner(pos, normal) + height;
    return dist;
}

float sdfTriPrism(const vec3 orgin, float h0, float h1) {
    //h[0] represents half the length of the base of the triangular prism
    //h[1] represents half the height of the prism along the z-axis
    vec3 pos;
    vec3_dup(pos, orgin);

    vec3 q = {fabsf(pos[0]), fabsf(pos[1]), fabsf(pos[2])};
    float dist = fmaxf(q[2] - h1, fmaxf(q[0] * 0.866025f + pos[1] * 0.5f, -pos[1]) - h0 * 0.5f);
    return dist;
}

float sdfVerticalCapsule(const vec3 origin, float height, float radius) {
    vec3 pos;
    vec3_dup(pos, origin);
    pos[1] -= clamp(pos[1], 0.0f, height);

    return vec3_len(pos) - radius;
}

float sdfSphere(const vec3 pos) {
    vec3 origin;
    vec3_set(origin, 0.0);
    vec3 displacement;
    vec3_sub(displacement, pos, origin);

    return vec3_len(displacement) - 3.5f;
}
