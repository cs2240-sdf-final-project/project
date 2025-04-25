#include <stdlib.h>
#include <stdio.h>
#include <stdio.h>
#include <assert.h>
#include <cctype>
#include <math.h>
#include "linmath.h"
#include <sys/stat.h>
#include <sys/types.h>

#include "hello.h"
#include "sim_random.h"
#include <iostream>
#include <vector>
#include <random>

int enzyme_dup;
int enzyme_dupnoneed;
int enzyme_out;
int enzyme_const;

// width of the scene band
const float distance_threshold = 5e-2f;

// we take steps at most this size in order to avoid missing
// sign changes in the directional derivatives
const float max_step_size = 1.0f;
// take this many steps to balance performance with exploring the entire scene
const int number_of_steps = 1'000;
// if our sdf gets smaller than this amount, we will consider it an intersection with the surface
const float contact_threshold = 1e-4f;

const int depth = 3;

const float lm_pi = 3.14159265358979323846f;

const float pathContinuationProb = 0.9f;

const int numberOfSampling = 10;

typedef struct{
    vec3 pos;
    vec3 dir;
    vec3 contribution;

}Segment;

enum BisectAffinity {
    BISECT_LEFT,
    BISECT_RIGHT,
    BISECT_STOP,
};

typedef BisectAffinity AFFINITY_FUNC(float evaluate_at, int iter_count, void *context);
typedef float MIDPOINT_FUNC(float a, float b, void *context);
inline float bisect(float t_min, float t_max,
    AFFINITY_FUNC *get_affinity, MIDPOINT_FUNC *get_midpoint, void *context) {
    int iter_count = 0;
    for (;;) {
        float t_mid = get_midpoint(t_min, t_max, context);
        BisectAffinity aff = get_affinity(t_mid, iter_count, context);
        if (aff == BISECT_STOP) {
            return t_mid;
        } else if (aff == BISECT_LEFT) {
            t_max = t_mid;
        } else if (aff == BISECT_RIGHT) {
            t_min = t_mid;
        }
        iter_count += 1;
    }
}

void ray_step(vec3 out, const vec3 origin, const vec3 direction, float t) {
    vec3 step_size;
    vec3_dup(step_size, direction);
    vec3_scale(step_size, step_size, t);
    vec3_add(out, origin, step_size);
}

float clamp(float x, float min, float max) {
    return fmaxf(fminf(x, max), min);
}

void dehomogenize(vec3 out, const vec4 in) {
    for (int i = 0; i < 3; i++) {
        out[i] = in[i] / in[3];
    }
}

float lerp(float x, float in_min, float in_max, float out_min, float out_max) {
    return out_min + (x - in_min) * (out_max - out_min) / (in_max - in_min);
}

void vec3_set(vec3 out, float value) {
    vec3 to_set = { value, value, value };
    vec3_dup(out, to_set);
}

void vec3_set(vec3 out, float v1, float v2, float v3) {
    vec3 to_set = { v1, v2, v3 };
    vec3_dup(out, to_set);
}

void vec2_abs(vec2 out, const vec2 in) {
    for (int i = 0; i < 2; i++) {
        out[i] = fabsf(in[i]);
    }
}

void vec3_cwiseProduct(vec3 out, vec3 a, vec3 b){
    out[0] = a[0]*b[0];
    out[1] = a[1]*b[1];
    out[2] = a[2]*b[2];
}

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

float sdfTriPrism(const vec3 origin, float h0, float h1) {
    //h[0] represents half the length of the base of the triangular prism
    //h[1] represents half the height of the prism along the z-axis
    vec3 pos;
    vec3_dup(pos, origin);

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

struct SceneParams {
    float object_1_x;
    float object_1_y;
    float object_1_z;
    float object_1_r;
    float object_1_h;
    float object_2_x;
    float object_2_y;
    float object_2_z;
    float object_2_r;
    float object_2_h;

    //light color
    float object_1_color[3];
    float object_1_direction[3];
    float object_color1;
};

int number_of_scene_params = (int)(sizeof(SceneParams) / sizeof(float));

SceneParams *make_scene_params() {
    SceneParams *out = (SceneParams *)calloc(sizeof(float), (size_t)number_of_scene_params);
    assert(out);
    return out;
}

void free_scene_params(SceneParams *params) {
    free(params);
}

static inline SceneParams *params_from_float_pointer(const float *params) {
    return (SceneParams *)params;
}

static const float *float_pointer_from_params(const SceneParams *out) {
    return (float *)out;
}

static float *float_pointer_from_params(SceneParams *out) {
    return (float *)out;
}

void scene_params_elementwise_add(SceneParams *out_params, const SceneParams *a, const SceneParams *b) {
    float *raw_out_params = float_pointer_from_params(out_params);
    const float *raw_a = float_pointer_from_params(a);
    const float *raw_b = float_pointer_from_params(b);
    for (int i = 0; i < number_of_scene_params; i++) {
        raw_out_params[i] = raw_a[i] + raw_b[i];
    }
}

void scene_params_elementwise_mul(SceneParams *out_params, const SceneParams *a, const SceneParams *b) {
    float *raw_out_params = float_pointer_from_params(out_params);
    const float *raw_a = float_pointer_from_params(a);
    const float *raw_b = float_pointer_from_params(b);
    for (int i = 0; i < number_of_scene_params; i++) {
        raw_out_params[i] = raw_a[i] * raw_b[i];
    }
}

void outer_product_add_assign(SceneParamsPerChannel *ppc, const SceneParams *params, const vec3 rgb) {
    const float *raw_params = float_pointer_from_params(params);
    for (int c = 0; c < 3; c++) {
        float *pc = float_pointer_from_params(ppc->rgb[c]);
        for (int i = 0; i < number_of_scene_params; i++) {
            pc[i] += raw_params[i] * rgb[c];
        }
    }
}

void scene_params_scale(SceneParams *out_params, const SceneParams *a, float scale_by) {
    float *raw_out_params = float_pointer_from_params(out_params);
    const float *raw_a = float_pointer_from_params(a);
    for (int i = 0; i < number_of_scene_params; i++) {
        raw_out_params[i] = raw_a[i] * scale_by;
    }
}
void scene_params_fill(SceneParams *params, float fill_with) {
    float *raw_out_params = float_pointer_from_params(params);
    for (int i = 0; i < number_of_scene_params; i++) {
        raw_out_params[i] = fill_with;
    }
}

float scene_parameter_get(const SceneParams *params, long p) {
    const float *raw_params = float_pointer_from_params(params);
    return raw_params[p];
}

void scene_params_set(SceneParams *params, long p, float value) {
    float *raw_params = float_pointer_from_params(params);
    raw_params[p] = value;
}

typedef struct {
    float distance;
    float ambient[3];
    float diffuse[3];
    float specular[3];
    float emissive[3];
    bool isReflected;
    float shininess;
} SdfResult;

static inline void default_scene_sample(SdfResult *s) {
    s->distance = INFINITY;
    vec3_set(s->ambient, 0.0);
    vec3_set(s->diffuse, 0.0);
    vec3_set(s->specular, 0.0);
    vec3_set(s->emissive, 0.0);
    s->isReflected = false;
    s->shininess = 1.0;
}

/** writes the output in the first paramemter */
static inline void compose_scene_sample(SdfResult *destination, SdfResult *b) {
    if (destination->distance < b->distance) {
        return; // nothing to be done
    }
    destination->distance = b->distance;
    vec3_dup(destination->ambient, b->ambient);
    vec3_dup(destination->diffuse, b->diffuse);
    vec3_dup(destination->specular, b->specular);
    vec3_dup(destination->emissive, b->emissive);
    destination->isReflected = b->isReflected;
    destination->shininess = b->shininess;
}

static inline void object_cylinder(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    vec3 doffset = {
        params->object_1_x,
        params->object_1_y,
        params->object_1_z};
    vec3 offset = {-0.4f, -0.2f, 0.2f};
    vec3_add(offset, offset, pos);
    vec3_add(offset, offset, doffset);
    sample->distance = sdfCylinder(offset, 0.3f + params->object_1_r, 0.5f + params->object_1_h);
    vec3_set(sample->diffuse, 0.4860f, 0.6310f, 0.6630f);
    vec3_set(sample->ambient, 0.4860f, 0.6310f, 0.6630f);
    vec3_set(sample->specular, 0.8f, 0.8f, 0.8f);
    //vec3_set(sample->emissive, 10.f, 10.f, 10.f);
    //sample->isReflected = true;
}
static inline void object_capsule(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    vec3 dpos = {
        params->object_2_x,
        params->object_2_y,
        params->object_2_z};
    vec3 offset = {0.4f, -0.3f, -0.5f};
    vec3_add(offset, offset, pos);
    vec3_add(offset, offset, dpos);
    sample->distance = sdfVerticalCapsule(offset, 0.5f, 0.3f);
    vec3_set(sample->diffuse, 0.4860f, 0.6310f, 0.6630f);
    vec3_set(sample->ambient, 0.4860f, 0.6310f, 0.6630f);
    //vec3_set(sample->emissive, 0.5f, 0.5f, 0.5f);
    vec3_set(sample->specular, 0.8f, 0.8f, 0.8f);
}
// Back:
static inline void object_backwall(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    (void)params;
    vec3 plane_normal = {0.0, 0.0, 1.0};
    sample->distance = sdfPlane(pos, plane_normal, 1.0f);
    vec3_set(sample->ambient, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->diffuse, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->specular, 0.4f);
}
// Ceiling:
static inline void object_topwall(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    (void)params;
    vec3 plane_normal = {0.0, 1.0, 0.0};
    sample->distance = sdfPlane(pos, plane_normal, 1.0f);
    vec3_set(sample->ambient, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->diffuse, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->specular, 0.4f);
    vec3_set(sample->emissive, 0.9f, 0.9f, 0.9f);
   
}
// Left:
static inline void object_leftwall(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    (void)params;
    vec3 plane_normal = {1.0, 0.0, 0.0};
    sample->distance = sdfPlane(pos, plane_normal, 1.0f);
    vec3_set(sample->ambient, 0.63f, 0.065f, 0.05f);
    vec3_set(sample->diffuse, 0.63f, 0.065f, 0.05f);
    vec3_set(sample->specular, 0.4f);
    sample->isReflected = true;
}
// Right:
static inline void object_rightwall(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    (void)params;
    vec3 plane_normal = {-1.0, 0.0, 0.0};
    sample->distance = sdfPlane(pos, plane_normal, 1.0f);
    vec3_set(sample->ambient, 0.14f, 0.45f, 0.091f);
    vec3_set(sample->diffuse, 0.14f, 0.45f, 0.091f);
    vec3_set(sample->specular, 0.4f);
}
// Floor:
static inline void object_bottomwall(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    (void)params;
    vec3 plane_normal = {0.0, -1.0, 0.0};
    sample->distance = sdfPlane(pos, plane_normal, 1.0f);
    vec3_set(sample->ambient, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->diffuse, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->specular, 0.4f);
    sample->isReflected = true;
}
inline void scene_sample(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    default_scene_sample(sample);
    SdfResult working;
    default_scene_sample(&working);
    object_cylinder(pos, params, &working);
    compose_scene_sample(sample, &working);
    default_scene_sample(&working);
    object_capsule(pos, params, &working);
    compose_scene_sample(sample, &working);
    default_scene_sample(&working);
    object_backwall(pos, params, &working);
    compose_scene_sample(sample, &working);
    default_scene_sample(&working);
    object_topwall(pos, params, &working);
    compose_scene_sample(sample, &working);
    default_scene_sample(&working);
    object_leftwall(pos, params, &working);
    compose_scene_sample(sample, &working);
    default_scene_sample(&working);
    object_rightwall(pos, params, &working);
    compose_scene_sample(sample, &working);
    default_scene_sample(&working);
    object_bottomwall(pos, params, &working);
    compose_scene_sample(sample, &working);
}

/** specialization of scene_sample that just returns the value of the sdf */
float scene_sample_sdf(const vec3 pos, const SceneParams *params) {
    SdfResult sample;
    scene_sample(pos, params, &sample);
    return sample.distance;
}

float directional_derivative_inner(const vec3 origin, const vec3 direction, float t, const SceneParams *params) {
    vec3 scaled;
    vec3_scale(scaled, direction, t);
    vec3 added;
    vec3_add(added, origin, scaled);
    return scene_sample_sdf(added, params);
}

extern float __enzyme_fwddiff_directional(void *,
    int, const float *,
    int, const float *,
    int, float, float,
    int, const SceneParams *);

float directional_derivative(const vec3 pos, const vec3 direction, float t, const SceneParams *params) {
    float dt = 1.0f;
    return __enzyme_fwddiff_directional(
        (void *)directional_derivative_inner,
        enzyme_const, pos,
        enzyme_const, direction,
        enzyme_dup, t, dt,
        enzyme_const, params);
}

typedef struct {
    vec3 origin, direction;
    const SceneParams *params;
    int iter_cap;
} CriticalPointContext;

static inline float critical_point_bisect_midpoint(float a, float b, void *context) {
    CriticalPointContext *ctx = (CriticalPointContext *)context;
    float dir_a = directional_derivative(ctx->origin, ctx->direction, a, ctx->params);
    float dir_b = directional_derivative(ctx->origin, ctx->direction, b, ctx->params);
    // have points (a, dir_a) and (b, dir_b). want to find the zero crossing.
    float denom = dir_a - dir_b;
    if (fabsf(denom) < 1e-5) {
        return a; // arbitrarily choose one of the endpoints
    }
    float numer = dir_a * b - a * dir_b;
    return numer / denom;
}

static inline BisectAffinity critical_point_bisect_affinity(float t, int iter_count, void *context) {
    CriticalPointContext *ctx = (CriticalPointContext *)context;
    if (iter_count >= ctx->iter_cap) {
        return BISECT_STOP;
    }
    // do a distance check on a fixed small iteration number to stop early
    const int early_stop_check_iter = 2;
    const float preliminary_distance_threshold = 1.0f;
    if (iter_count == early_stop_check_iter) {
        vec3 pos;
        ray_step(pos, ctx->origin, ctx->direction, t);
        if (scene_sample_sdf(pos, ctx->params) > preliminary_distance_threshold) {
            return BISECT_STOP;
        }
    }
    float dir_t = directional_derivative(ctx->origin, ctx->direction, t, ctx->params);
    bool is_approaching = dir_t < 0;
    return is_approaching ? BISECT_RIGHT : BISECT_LEFT;
}

typedef struct {
    bool found_critical_point;
    float t_if_found_critical_point;
} SearchResult;

void search_for_critical_point(const vec3 origin, const vec3 direction, const SceneParams *params, float t_min, float t_max, SearchResult *ret) {
    ret->found_critical_point = false;
    ret->t_if_found_critical_point = 0.0;

    // we consider directional derivatives this large to be zero
    const float directional_derivative_threshold = 1e-1f;
    // we will bisect this many iterations in order to find the true location of the minimum directional derivative
    const int extensive_depth = 12;
    CriticalPointContext context;
    vec3_dup(context.origin, origin);
    vec3_dup(context.direction, direction);
    context.params = params;
    context.iter_cap = extensive_depth;
    float best_t = bisect(t_min, t_max, critical_point_bisect_affinity, critical_point_bisect_midpoint, &context);
    vec3 best_pos;
    ray_step(best_pos, origin, direction, best_t);
    if (scene_sample_sdf(best_pos, params) > distance_threshold) {
        return;
    }
    float best_dir_t = directional_derivative(origin, direction, best_t, params);
    if (fabsf(best_dir_t) > directional_derivative_threshold) {
        return;
    }
    ret->found_critical_point = true;
    ret->t_if_found_critical_point = best_t;
}

float sdf_normal_wrapper(const vec3 pos, const float *params) {
    SceneParams* scene_params = params_from_float_pointer(params);
    return scene_sample_sdf(pos, scene_params);
}

extern void __enzyme_autodiff_normal(void *, int, const float *, float *, int, const float *);

void get_normal_from(vec3 normal, const vec3 pos, const SceneParams *params) {
    vec3 dpos;
    vec3_set(dpos, 0.0f);
    const float *raw_params = float_pointer_from_params(params);
    __enzyme_autodiff_normal(
        (void*)sdf_normal_wrapper,
        enzyme_dup, pos, dpos,
        enzyme_const, raw_params
    );
    vec3_dup(normal, dpos);
}

float sdf_theta_wrapper(const vec3 pos, const float *params) {
    SceneParams *scene_params = params_from_float_pointer(params);
    return scene_sample_sdf(pos, scene_params);
}

extern void __enzyme_autodiff_theta(void *, int,const float *, int, const float *, float *);

void diff_sdf(const vec3 pos, SceneParams *paramsOut, const SceneParams *paramsIn) {
    const float *raw_params = float_pointer_from_params(paramsIn);
    float* draw_params = float_pointer_from_params(paramsOut);
    __enzyme_autodiff_theta(
        (void*)sdf_theta_wrapper,
        enzyme_const, pos,
        enzyme_dup, raw_params, draw_params
    );
}

void phongLight(vec3 radiance, const vec3 looking, const vec3 normal, const SdfResult *sample, const SceneParams *params) {
    float lightColors[3][3] = {
        {.8f, .8f, .8f},
        {.2f, .2f, .2f},
        {.2f, .2f, .2f},
    };

    float lightDirections[3][3] = {
        {0.f, -1.f, 0.f},
        {-3.f, 2.f, 0.f},
        {0.f, -1.f, 3.f},
    };
    float kd = 0.5;
    float ks = 1.0;
    vec3_dup(radiance, sample->ambient);
    for (int l = 0; l < 3; l++) {
        vec3 lightColor; 
        vec3_dup(lightColor, lightColors[l]);
        lightColor[0] += params->object_1_color[0];
        lightColor[1] += params->object_1_color[1]; 
        lightColor[2] += params->object_1_color[2]; 
        vec3 light_dir;
        vec3_norm(light_dir, lightDirections[l]);
        light_dir[0] += params->object_1_direction[0];
        light_dir[1] += params->object_1_direction[1];
        light_dir[2] += params->object_1_direction[2];
        float facing = fmaxf(0.0, vec3_mul_inner(normal, light_dir));

        
        vec3 bounce;
        for (int i = 0; i < 3; i++) {
            bounce[i] = light_dir[i] - normal[i] * facing * 2.f;
        }
        float specular = powf(vec3_mul_inner(bounce, looking), sample->shininess);
        for (int i = 0; i < 3; i++) {
            radiance[i] += kd * facing * sample->diffuse[i] * lightColor[i];
            radiance[i] += ks * specular * sample->specular[i] * lightColor[i];
        }
    }
}

typedef struct {
    bool found_intersection;
    float intersection_t;
} IntersectionResult;

IntersectionResult trace_ray(

    const vec3 origin,
    const vec3 direction,
    const SceneParams *params

){
    float t = 0.0;

    IntersectionResult ret = { false, 0.0 };
    for (int i = 0; i < number_of_steps; i++) {

        vec3 pos;
        ray_step(pos, origin, direction, t);
        float distance = scene_sample_sdf(pos, params);
        if (distance < contact_threshold) {
            ret.found_intersection = true;
            ret.intersection_t = t;
            break;
        }
        float step_size = fminf(max_step_size, distance);
        t += step_size;
    }
    return ret;

}

IntersectionResult trace_ray_get_critical_point(
    SearchResult *critical_point,
    const vec3 origin,
    const vec3 direction,
    const SceneParams *params
) {


    float t = 0.0;
    float previous_t = 0.0;
    critical_point->found_critical_point = false;
    IntersectionResult ret = { false, 0.0 };
    for (int i = 0; i < number_of_steps; i++) {
        float dir_previous_t = directional_derivative(origin, direction, previous_t, params);
        float dir_t = directional_derivative(origin, direction, t, params);
        // look for a sign change in the directional derivative
        if (!critical_point->found_critical_point && ((dir_previous_t < 0) && (dir_t > 0))) {
            // let's try to find critical_point between t and previous_t;
            search_for_critical_point(origin, direction, params, previous_t, t, critical_point);
        }
        vec3 pos;
        ray_step(pos, origin, direction, t);
        float distance = scene_sample_sdf(pos, params);
        if (distance < contact_threshold) {
            ret.found_intersection = true;
            ret.intersection_t = t;
            break;
        }
        float step_size = fminf(max_step_size, distance);
        previous_t = t;
        t += step_size;
    }
    return ret;
}

void uniformSampleHemisphere(vec3 direction, vec3 normal){

    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(0.0f, 1.0f);

    float r1 = dis(gen);
    float r2 = dis(gen);

    float phi = 2.0f * lm_pi * r1;
    float cosTheta = 1.0f - r2;
    //float cosTheta = r2;
    float sinTheta = std::sqrt(1.0f - cosTheta * cosTheta);

    float x = sinTheta * std::cos(phi);
    float y = sinTheta * std::sin(phi);
    float z = cosTheta;

    vec3 wi;
    vec3_dup(wi,vec3{x,y,z});

    vec3 v1;
    vec3_dup(v1,vec3{0.0f, 0.0f, 1.0f});
    vec3 v2;
    vec3_dup(v2,normal);
    vec3 v;
    vec3_dup(v,wi);
    vec3 n_normal;
    vec3_norm(n_normal, v2);

    vec3 k;
    vec3_mul_cross(k,v1,n_normal);

    float dot_product = vec3_mul_inner(v1,n_normal);

    vec3 wi_rotated;
    vec3_set(wi_rotated, 0.0f);

    float epsilon = 1e-6f;

    if(dot_product > 1.0f - epsilon){

        vec3_dup(direction, wi); 
        return;

    }else if(dot_product < -1.0f + epsilon){

        vec3_dup(wi_rotated , vec3{wi[0], -wi[1], -wi[2]});
        vec3_dup(direction, wi_rotated);

        return;
    }

    vec3 k_norm;
    vec3_norm(k_norm,k);

    float costheta = vec3_mul_inner(v1,n_normal);

    costheta = fmaxf(-1.0f, fminf(1.0f, costheta));
  
    float sintheta = vec3_len(k);

    
    vec3 first_term;
    vec3_scale(first_term, wi, costheta);
    vec3 second_term;
    vec3_mul_cross(second_term,k_norm,wi);
    vec3_scale(second_term, second_term,sintheta);

    vec3 last_term;
    float dotk_w = vec3_mul_inner(k_norm,wi);
    dotk_w *= (1 - costheta);
    vec3_scale(last_term, k_norm,dotk_w);

    vec3_add(wi_rotated, wi_rotated,first_term);
    vec3_add(wi_rotated, wi_rotated,second_term);
    vec3_add(wi_rotated, wi_rotated,last_term);





    vec3_dup(direction,wi_rotated);





}
void test_hemisphere_sampling() {
    printf("Testing uniformSampleHemisphere orientation...\n");
    // Adjust num_samples as needed
    int num_samples = 10000; // How many random directions to test per case
    int error_count = 0;
    // Use a small negative tolerance for checks like >= 0
    // to account for minor floating point inaccuracies near zero
    float tolerance = -1e-6f;

    // --- Test Case 1: Normal = (0, 0, 1) ---
    printf(" Test 1: Normal = (0, 0, 1). Expecting dir[2] >= 0\n");
    vec3 test_normal_z_pos;
    vec3_dup(test_normal_z_pos, vec3{0.f, 0.f, 1.f}); // Use the test normal
    for (int i = 0; i < num_samples; ++i) {
        vec3 test_dir;
        // --- Call your function correctly ---
        // If it needs PDF: float pdf; uniformSampleHemisphere(test_dir, pdf, test_normal_z_pos);
        uniformSampleHemisphere(test_dir, test_normal_z_pos); // Pass the correct normal
        // --- Check the result ---
        if (test_dir[2] < tolerance) { // Check z component >= 0 (allowing for tolerance)
            printf("  ERROR (Sample %d): Normal=(0,0,1), Dir=(%f, %f, %f) -> dir[2] < 0!\n",
                   i, test_dir[0], test_dir[1], test_dir[2]);
            error_count++;
             // break; // Optional: Stop this test case on first error
        }
    }

    // --- Test Case 2: Normal = (0, 1, 0) ---
    printf(" Test 2: Normal = (0, 1, 0). Expecting dir[1] >= 0\n");
    vec3 test_normal_y_pos;
    vec3_dup(test_normal_y_pos, vec3{0.f, 1.f, 0.f}); // Use the test normal
    for (int i = 0; i < num_samples; ++i) {
        vec3 test_dir;
        uniformSampleHemisphere(test_dir, test_normal_y_pos); // Pass the correct normal
        if (test_dir[1] < tolerance) { // Check y component >= 0
            printf("  ERROR (Sample %d): Normal=(0,1,0), Dir=(%f, %f, %f) -> dir[1] < 0!\n",
                   i, test_dir[0], test_dir[1], test_dir[2]);
            error_count++;
            // break;
        }
    }

    // --- Test Case 3: Normal = (1, 0, 0) ---
    printf(" Test 3: Normal = (1, 0, 0). Expecting dir[0] >= 0\n");
    vec3 test_normal_x_pos;
    vec3_dup(test_normal_x_pos, vec3{1.f, 0.f, 0.f}); // Use the test normal
    for (int i = 0; i < num_samples; ++i) {
        vec3 test_dir;
        uniformSampleHemisphere(test_dir, test_normal_x_pos); // Pass the correct normal
        if (test_dir[0] < tolerance) { // Check x component >= 0
            printf("  ERROR (Sample %d): Normal=(1,0,0), Dir=(%f, %f, %f) -> dir[0] < 0!\n",
                   i, test_dir[0], test_dir[1], test_dir[2]);
            error_count++;
            // break;
        }
    }

    // --- Test Case 4: Normal = (0, 0, -1) ---
    printf(" Test 4: Normal = (0, 0, -1). Expecting dir[2] <= 0\n");
    vec3 test_normal_z_neg;
    vec3_dup(test_normal_z_neg, vec3{0.f, 0.f, -1.f}); // Use the test normal
    for (int i = 0; i < num_samples; ++i) {
        vec3 test_dir;
        uniformSampleHemisphere(test_dir, test_normal_z_neg); // Pass the correct normal
        // Check z component <= 0 (allowing for tolerance)
        // test_dir[2] > -tolerance is equivalent to test_dir[2] > tolerance_positive
        if (test_dir[2] > -tolerance) {
            printf("  ERROR (Sample %d): Normal=(0,0,-1), Dir=(%f, %f, %f) -> dir[2] > 0!\n",
                   i, test_dir[0], test_dir[1], test_dir[2]);
            error_count++;
            // break;
        }
    }

    // --- Summary ---
    printf("----------------------------------------\n");
    if (error_count == 0) {
        printf("Hemisphere sampling orientation test PASSED for %d samples per case.\n", num_samples);
    } else {
        printf("Hemisphere sampling orientation test FAILED with %d error(s).\n", error_count);
    }
    printf("----------------------------------------\n");
}

std::vector<Segment> getSecondaryPath(vec3 origin, vec3 direction,RandomState *random,const SceneParams *params){

    std::vector<Segment> path;

    vec3 current_position;
    vec3_dup(current_position, origin);
    vec3 current_direction;
    vec3_dup(current_direction,direction);
    // IntersectionResult hitPoint = trace_ray(current_position,current_direction,params);
    
    // //store current ray 
    // Segment newSegment;
    // vec3_dup(newSegment.pos,current_position);
    // vec3_dup(newSegment.dir,current_direction);
    // path.push_back(newSegment);
    

    while(true){ 

        IntersectionResult hitPoint = trace_ray(current_position,current_direction,params);
        if(!hitPoint.found_intersection){
            return path;
        }
        //store current ray 
        Segment newSegment;
        vec3_dup(newSegment.pos,current_position);
        vec3_dup(newSegment.dir,current_direction);
        path.push_back(newSegment);

        static thread_local std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<float> dis(0.0f, 1.0f);

        if(random_next_float(random)> pathContinuationProb){
            break;
        }

        //find intersection point

        vec3 hit_intersection_point;
        ray_step(hit_intersection_point, current_position, current_direction, hitPoint.intersection_t);

        SdfResult sample;
        scene_sample(hit_intersection_point, params, &sample);

        vec3 normal;
        get_normal_from(normal, hit_intersection_point, params);


        vec3 newDir;
        if(sample.isReflected){
            vec3 refl_dir;
            vec3_reflect(refl_dir,current_direction,normal);
            vec3_dup(newDir,refl_dir);
            //sample_hemisphere(newDir, normal, random);
        }else{

            //uniformSampleHemisphere(newDir, normal);
            sample_hemisphere(newDir, normal, random);

        }
        //sample_hemisphere(newDir, normal, random);
        //uniformSampleHemisphere(newDir, normal);

        const float ray_epsilon = 1e-4f;


        vec3 offset_position;
        vec3_scale(offset_position, normal, ray_epsilon);
        vec3_add(hit_intersection_point,hit_intersection_point,offset_position);


        //update new ray 
        vec3_dup(current_position,hit_intersection_point);
        vec3_dup(current_direction,newDir);


    }

    return path;


    
}


void get_radiance_at(
    vec3 radiance,
    std::vector<Segment> &path,
    const SceneParams *params,
    const int depth,
    RandomState* random

) {
    vec3 spectralFilter;
    vec3_set(spectralFilter, 1.f);
    vec3 intensity;
    vec3_set(intensity, 0.f);
    // printf("Path Start Position: (%f, %f, %f)\n",
    //            path[0].pos[0],
    //            path[0].pos[1],
    //            path[0].pos[2]);


    
    for(size_t i = 0; i < path.size(); i++){

        if(i != path.size() - 1){

            vec3 hitPosition;
            vec3_dup(hitPosition, path[i+1].pos);
            vec3 hitDirection;
            vec3_dup(hitDirection,path[i].dir);
            vec3 wi;
            vec3_dup(wi,path[i+1].dir);

            

            SdfResult sample;
            scene_sample(hitPosition, params, &sample);
            vec3 normal;
            get_normal_from(normal, hitPosition, params);

            // float normal_flip = vec3_mul_inner(normal,hitDirection);

            // if(normal_flip > 0){
            //     vec3_scale(normal,normal,-1.f);
            // }

           
            // vec3 normal_color;
            // vec3_scale(normal_color, normal, 0.5f);
            // vec3_add(normal_color, normal_color, vec3{0.5f, 0.5f, 0.5f});
            // vec3_dup(radiance, normal_color);
            
                        

            vec3 emissive;
            vec3_dup(emissive, sample.emissive);


            vec3 emissive_part;
            vec3_cwiseProduct(emissive_part, emissive, spectralFilter);
           
      
            
            vec3_add(intensity, intensity, emissive_part);
            //vec3_add(intensity,intensity,vec3{0.1f,0.1f,0.1f});

            if(sample.isReflected){
                vec3_scale(spectralFilter,spectralFilter, 1.0f / pathContinuationProb);
                vec3 brdf;
                vec3_set(brdf, 1.0f);
                vec3_cwiseProduct(spectralFilter, spectralFilter, brdf);


            }else{
                vec3_scale(spectralFilter,spectralFilter, 1.0f / pathContinuationProb);
      

                float pdf = 2.f * lm_pi;
                vec3_scale(spectralFilter,spectralFilter, pdf);

                float cosine_term = fmaxf(0.f, vec3_mul_inner(normal, wi));
                    //float cosine_term = vec3_mul_inner(normal, wi);
                vec3_scale(spectralFilter, spectralFilter, cosine_term);

                vec3 brdf;
                vec3_scale(brdf, sample.diffuse, 1.0f / lm_pi);
                vec3_cwiseProduct(spectralFilter, spectralFilter, brdf);
                

            }
             

           

        }
        
    }

    vec3_dup(radiance,intensity);
    // printf("  i=%zu: intensity (Pre-Emissive): (%f, %f, %f)\n",
    //                 intensity[0], intensity[1], intensity[2]);
   
    
    
    
}


GradientImage make_gradient_image(long image_width, long image_height) {
    const long num_subpixels = 3;
    long num_floats = number_of_scene_params * image_width * image_height * num_subpixels;
    float *buf = (float *)calloc(sizeof(float), (size_t)num_floats);
    assert(buf);
    GradientStrides strides;
    strides.parameter_stride = 1;
    strides.subpixel_stride = number_of_scene_params;
    strides.col_stride = number_of_scene_params * num_subpixels;
    strides.row_stride = number_of_scene_params * num_subpixels * image_width;
    GradientImage ret;
    ret.strides = strides;
    ret.image_height = image_height;
    ret.image_width = image_width;
    ret.num_subpixels = 3;
    ret.buf = buf;
    return ret;
}

void gradient_image_slice(Image *image, const GradientImage *gradient, long parameter_no) {
    SceneParamsPerChannel ppc;
    for (int ch = 0; ch < 3; ch++) {
        ppc.rgb[ch] = make_scene_params();
    }
    for (long r = 0; r < image->image_height; r++) {
        for (long c = 0; c < image->image_height; c++) {
            gradient_image_get(&ppc, gradient, r, c);

            vec3 radiance;
            for (int ch = 0; ch < 3; ch++) {
                float param_value = scene_parameter_get(ppc.rgb[ch], parameter_no);
                radiance[ch] = param_value;
            }
            vec3 half;
            vec3_set(half, 0.5f);
            vec3_scale(radiance, radiance, 0.5f);
            vec3_add(radiance, radiance, half);
            image_set(image, r, c, radiance);
        }
    }
    for (int ch = 0; ch < 3; ch++) {
        free_scene_params(ppc.rgb[ch]);
    }
}

void free_gradient_image(GradientImage *image) {
    free(image->buf);
}

static inline long gradient_image_get_index(const GradientStrides *s, long r, long c, long subpixel, long param) {
    return r * s->row_stride + c * s->col_stride + subpixel * s->subpixel_stride + param * s->parameter_stride;
}

void gradient_image_set(const SceneParamsPerChannel *ppc, GradientImage *image, long ir, long ic) {
    for (long subpixel = 0; subpixel < 3; subpixel++) {
        const SceneParams *params = ppc->rgb[subpixel];
        const float *raw_params = float_pointer_from_params(params);
        for (long p = 0; p < number_of_scene_params; p++) {
            long index = gradient_image_get_index(&image->strides, ir, ic, subpixel, p);
            image->buf[index] = raw_params[p];
        }
    }
}

void gradient_image_get(SceneParamsPerChannel *ppc, const GradientImage *image, long ir, long ic) {
    for (long subpixel = 0; subpixel < 3; subpixel++) {
        SceneParams *params = ppc->rgb[subpixel];
        float *raw_params = float_pointer_from_params(params);
        for (long p = 0; p < number_of_scene_params; p++) {
            long index = gradient_image_get_index(&image->strides, ir, ic, subpixel, p);
            raw_params[p] = image->buf[index];
        }
    }
}

void render_get_radiance_wrapper(
    vec3 radiance,
    std::vector<Segment> &path,
    const float *raw_params,
    RandomState* random
) {
    const SceneParams *params = params_from_float_pointer(raw_params);

    get_radiance_at(radiance, path, params,depth,random);
}

extern void __enzyme_fwddiff_radiance(
    void *,
    int, float *, float *,
    int, std::vector<Segment>&,
    int, const float *, const float *,
    int, RandomState *
);

void render_pixel(
    vec3 real,
    const vec3 origin,
    const vec3 direction,
    SceneParamsPerChannel *params_per_channel,
    const SceneParams *params,
    RandomState* random
) {
    SceneParams *dummy_params = make_scene_params();
    SearchResult critical_point;
    IntersectionResult intersection = trace_ray_get_critical_point(&critical_point, origin, direction, params);

    
    vec3 o;
    vec3_dup(o,origin);
    vec3 d;
    vec3_dup(d,direction);

    std::vector<Segment> path = getSecondaryPath(o , d ,random,params);

    for (int p = 0; p < number_of_scene_params; p++) {
        vec3 radiance;
        vec3_set(radiance, 1.f);
        vec3 d_radiance;
        vec3_set(d_radiance, 1.f);

        const float *raw_params = float_pointer_from_params(params);
        scene_params_fill(dummy_params, 0.f);
        scene_params_set(dummy_params, p, 1.f);
        const float *raw_dummy_params = float_pointer_from_params(dummy_params);

        __enzyme_fwddiff_radiance(
            (void*)render_get_radiance_wrapper,
            enzyme_dup, radiance, d_radiance,
            enzyme_const, path,
            enzyme_dup, raw_params, raw_dummy_params,
            enzyme_const, random
        );

        for (int ch = 0; ch < 3; ch++) {
            scene_params_set(params_per_channel->rgb[ch], p, d_radiance[ch]);
        }
    }

    get_radiance_at(real, path, params,depth,random);

    // if(critical_point.found_critical_point) {
    //     vec3 y_star;
    //     ray_step(y_star, origin, direction, critical_point.t_if_found_critical_point);
    //     SdfResult sample;
    //     scene_sample(y_star, params, &sample);
    //     vec3 normal;
    //     get_normal_from(normal, y_star, params);
    //     vec3 y_star_radiance;
    //     phongLight(y_star_radiance, direction, normal, &sample,params);
    //     diff_sdf(y_star, dummy_params, params);
    //     vec3 deltaL;
    //     vec3_sub(deltaL, y_star_radiance, real);
    //     vec3_scale(deltaL, deltaL, 1 / distance_threshold);
    //     outer_product_add_assign(params_per_channel, dummy_params, deltaL);
    // }
    free_scene_params(dummy_params);
}

void render_pixel_wrapper(
    vec3 real,
    const vec3 origin,
    const vec3 direction,
    SceneParamsPerChannel *params_per_channel,
    const SceneParams *params,
    RandomState* random
) { 
    
    for(int i = 0; i < numberOfSampling; i++){
        vec3 temp_real;
        render_pixel(temp_real, origin,direction,params_per_channel,params,random);
        vec3_add(real,real, temp_real);
    }
    vec3_scale(real,real, 1.f/numberOfSampling);

}



long get_index(Strides *s, long r, long c, long p) {
    return r * s->row_stride + c * s->col_stride + p * s->subpixel_stride;
}

//////////////////////////////////////
// Begin PPM Parser from ChatGPT /////
void skip_whitespace_and_comments(FILE *f) {
    int c;
    while ((c = fgetc(f)) != EOF) {
        if (c == '#') {
            // Skip the comment line
            while ((c = fgetc(f)) != '\n' && c != EOF);
        } else if (!isspace(c)) {
            ungetc(c, f);
            break;
        }
    }
}
int image_read_bpm(Image *image, FILE *f) {
    char header[3];
    if (fscanf(f, "%2s", header) != 1 || strcmp(header, "P6") != 0) {
        fprintf(stderr, "Unsupported or invalid PPM format\n");
        return 0;
    }
    skip_whitespace_and_comments(f);
    if (fscanf(f, "%ld", &image->image_width) != 1) return 0;

    skip_whitespace_and_comments(f);
    if (fscanf(f, "%ld", &image->image_height) != 1) return 0;

    skip_whitespace_and_comments(f);
    int maxval;
    if (fscanf(f, "%d", &maxval) != 1 || maxval != 255) {
        fprintf(stderr, "Unsupported max color value (only 255 supported)\n");
        return 0;
    }
    // Skip the single whitespace character after maxval
    fgetc(f);
    image->num_bytes = image->image_width * image->image_height * 3;
    image->buf = (char *)malloc((size_t)image->num_bytes);
    if (!image->buf) {
        fprintf(stderr, "Memory allocation failed\n");
        return 0;
    }
    size_t bytes_read = fread(image->buf, 1, (size_t)image->num_bytes, f);
    if (bytes_read != (size_t)image->num_bytes) {
        fprintf(stderr, "Failed to read image data\n");
        free(image->buf);
        image->buf = NULL;
        return 0;
    }
    return 1; // success
}
// End PPM Parser from ChatGPT ///////
//////////////////////////////////////

Image make_image(long image_width, long image_height) {
    Image image;
    long num_bytes = image_width * image_height * 3;
    assert(num_bytes > 0);
    char *buf = (char*)malloc((size_t)num_bytes);
    assert(buf);
    Strides strides;
    strides.row_stride = image_width * 3;
    strides.col_stride = 3;
    strides.subpixel_stride = 1;
    image.strides = strides;
    image.image_width = image_width;
    image.image_height = image_height;
    image.buf = buf;
    image.num_bytes = num_bytes;
    return image;
}

void free_image(Image *image) {
    free(image->buf);
}

// https://nullprogram.com/blog/2017/11/03/
void image_write_ppm(Image *image, FILE *f) {
    fprintf(f, "P6\n%ld %ld\n255\n", image->image_width, image->image_height);
    assert(image->num_bytes > 0);
    fwrite(image->buf, 1, (size_t)image->num_bytes, f);
    fflush(f);
}

void image_set(Image *image, long ir, long ic, const vec3 radiance) {
    assert(ir >= 0);
    assert(ir < image->image_height);
    assert(ic >= 0);
    assert(ic < image->image_width);
    for (long p = 0; p < 3; p++) {
        long index = get_index(&image->strides, ir, ic, p);
        char value = (char)clamp(radiance[p] * 255.f, 0.f, 255.f);
        image->buf[index] = value;
    }
}

void image_get(vec3 radiance, Image *image, long ir, long ic) {
    assert(ir >= 0);
    assert(ir < image->image_height);
    assert(ic >= 0);
    assert(ic < image->image_width);
    for (long p = 0; p < 3; p++) {
        long index = get_index(&image->strides, ir, ic, p);
        char value = image->buf[index];
        radiance[p] = (float)value / 255.f;
    }
}

void render_image(Image *real, GradientImage *gradient, const SceneParams *params,RandomState *random) {
    float aspect = (float)real->image_width / (float)real->image_height ;
    float near_clip = 0.1f;
    float far_clip = 100.0f;
    //float y_fov = 0.785f;
    float y_fov = 1.f;
    vec3 camera_position = {0.f, 0.f, 3.6f};
    vec3 center = {0};
    vec3 up = {0.0, 1.0, 0.0};

    mat4x4 look_at;
    mat4x4_look_at(look_at, camera_position, center, up);

    mat4x4 perspective;
    mat4x4_perspective(perspective, y_fov, aspect, near_clip, far_clip);

    mat4x4 projection_view;
    mat4x4_mul(projection_view, perspective, look_at);

    mat4x4 world_from_camera;
    mat4x4_invert(world_from_camera, projection_view);

    // SceneParamsPerChannel ppc;
    // for (int c = 0; c < 3; c++) {
    //     ppc.rgb[c] = make_scene_params();
    // }
    #pragma omp parallel for 

    for (long ir = 0; ir < real->image_height; ir++) {
        // if (ir % 50 == 0) {
        //     printf("on row %ld\n", ir);
        // }
        SceneParamsPerChannel ppc;
        for (int c = 0; c < 3; c++) {
            ppc.rgb[c] = make_scene_params();
        }
        for (long ic = 0; ic < real->image_width; ic++) {
            float r = (float)ir;
            float c = (float)ic;

            float device_x = lerp(c, 0.0, (float)real->image_width, -1.0, 1.0);
            float device_y = lerp(r, 0.0, (float)real->image_height, -1.0, 1.0);

            vec4 unprojected = {device_x, device_y, 1.0, 1.0};
            vec4 homo;
            mat4x4_mul_vec4(homo, world_from_camera, unprojected);
            vec3 far;
            dehomogenize(far, homo);

            vec3 direction;
            vec3_sub(direction, far, camera_position);
            vec3_norm(direction, direction);

            // Calculate radiance and gradients for a single pixel
            vec3 out_real;
            //render_pixel(out_real, camera_position, direction, &ppc, params,random);
            render_pixel_wrapper(out_real, camera_position, direction, &ppc, params,random);

            image_set(real, ir, ic, out_real);
            gradient_image_set(&ppc, gradient, ir, ic);
        }
        for (int c = 0; c < 3; c++) {
        free_scene_params(ppc.rgb[c]);
        }
    }

    // for (int c = 0; c < 3; c++) {
    //     free_scene_params(ppc.rgb[c]);
    // }

    //test_hemisphere_sampling();
}
