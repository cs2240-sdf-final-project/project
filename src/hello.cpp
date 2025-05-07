#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <stdio.h>
#include <assert.h>
#include <cctype>
#include <math.h>
#include "linmath.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <inttypes.h>
#include <functional>

#include "hello.h"
#include "sim_random.h"
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

int enzyme_dup;
int enzyme_dupnoneed;
int enzyme_out;
int enzyme_const;

// blur kernel effect radius
const int blur_radius = 5;

// we take steps at most this size in order to avoid missing
// sign changes in the directional derivatives
const float max_step_size = 1.0f;
// take this many steps to balance performance with exploring the entire scene
const int number_of_steps = 1'00;
// if our sdf gets smaller than this amount, we will consider it an intersection with the surface
const float contact_threshold = 1e-4f;
// width of the scene band
const float distance_threshold = 1e-2f;

const float lm_pi = 3.14159265358979323846f;

const float pathContinuationProb = 0.9f;
// finite differences half epsilon
const float finite_difference_epsilon = 1e-3f;

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
    vec3_scale(step_size, direction, t);
    vec3_add(out, origin, step_size);
}

long clampl(long x, long min, long max) {
    return std::max(std::min(x, max), min);
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

float vec3_distance(const vec3 a, const vec3 b) {
    vec3 displacement;
    vec3_sub(displacement, b, a);
    return vec3_len(displacement);
}

void vec3_set(vec3 out, float value) {
    vec3 to_set = { value, value, value };
    vec3_dup(out, to_set);
}

void vec3_componentwise_mul(vec3 out, const vec3 a, const vec3 b) {
    for (int i = 0; i < 3; i++) {
        out[i] = a[i] * b[i];
    }
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

long lclamp(long x, long min, long max) {
    if (x > max) {
        return max;
    }
    if (x < min) {
        return min;
    }
    return x;
}

// maps 0..=dim-1 points to -1.0..=1.0
inline void translate_grid(vec3 out, const long point[3], const long dim[3]) {
    for (int i = 0; i < 3; i++) {
        out[i] = lerp((float)point[i], 0.0, (float)(dim[i] -1), -1.0, 1.0);
    }
}

// maps 0..=dim-1 points to -1.0..=1.0
inline void inv_translate_grid(float out[3], const float point[3], const long dim[3]) {
    for (int i = 0; i < 3; i++) {
        out[i] = lerp(point[i], -1.0, 1.0, 0.0, (float)(dim[i] -1));
    }
}

inline void sample_existing_sdf(
    const long strides[3],
    const long dim[3],
    float *sds,
    std::function<float(const vec3 pos)> sampler
) {
    for (long x0 = 0; x0 < dim[0]; x0++) {
        for (long x1 = 0; x1 < dim[1]; x1++) {
            for (long x2 = 0; x2 < dim[2]; x2++) {
                long x[3] = {x0, x1, x2};
                vec3 fx;
                translate_grid(fx, x, dim);
                long i = x0 * strides[0] + x1 * strides[1] + x2 * strides[2];
                assert(i >= 0 && i < dim[0] * dim[1] * dim[2]);
                sds[i] = sampler(fx);
            }
        }
    }
}

inline float sdfGrid(
    const vec3 pos,
    const long strides[3],
    const long dim[3],
    const float *sds
) {
    // calculate floor grid point
    float pos_in_grid[3];
    inv_translate_grid(pos_in_grid, pos, dim);
    long top_left_grid[3];
    for (long d = 0; d < 3; d++) {
        top_left_grid[d] = static_cast<long>(pos_in_grid[d]);
    }
    float weighted_sum = 0.0;
    float total_weights = 0.0;
    for (long dix = 0; dix < 2; dix++) {
        for (long diy = 0; diy < 2; diy++) {
            for (long diz = 0; diz < 2; diz++) {
                long di[3] = {dix, diy, diz};
                float weight = 1.0f;
                for (long d = 0; d < 3; d++) {
                    float w = pos_in_grid[d] - (float)top_left_grid[d];
                    weight *= di[d] == 0 ? 1.0f - w : w;
                }
                long index = 0;
                for (long d = 0; d < 3; d++) {
                    long coord = lclamp(top_left_grid[d] + di[d], 0, dim[d] - 1);
                    index += coord * strides[d];
                }
                assert(index >= 0 && index < dim[0] * dim[1] * dim[2]);
                float distance = sds[index];
                weighted_sum += distance * weight;
                total_weights += weight;
            }
        }
    }
    return weighted_sum / total_weights;
}

inline void normalGrid(
    vec3 normal_out,
    const long pos[3],
    const long strides[3],
    const long dim[3],
    const float *sds
) {
    float step_size[3];
    for (int axis = 0; axis < 3; axis++) {
        step_size[axis] = 1.0f / static_cast<float>(dim[axis] - 1);
    }
    long center_index = 0;
    for (long d = 0; d < 3; d++) {
        assert(pos[d] >= 0 && pos[d] < dim[d]);
        center_index += pos[d] * strides[d];
    }
    float center_value = sds[center_index];
    for (long axis = 0; axis < 3; axis++) {
        float normal_sum = 0.0;
        float normal_count = 0.0;
        // visit dd = -1, dd = 1
        for (long dd = -1; dd <= 1; dd += 2) {
            if (pos[axis] + dd < 0 || pos[axis] + dd >= dim[axis]) {
                continue;
            }
            long sample_index = center_index + dd * strides[axis];
            float sample_value = sds[sample_index];
            float diff = sample_value - center_value;
            normal_sum += diff * static_cast<float>(dd);
            normal_count += step_size[axis];
        }
        normal_out[axis] = normal_sum / normal_count;
        // printf("%g\n", normal_out[axis]);
    }
}

inline float grid_consistency_loss(
    const long strides[3],
    const long dim[3],
    const float *sds
) {
    float total = 0.f;
    for (long x = 0; x < dim[0]; x++) {
        for (long y = 0; y < dim[1]; y++) {
            for (long z = 0; z < dim[2]; z++) {
                long xyz[3] = {x, y, z};
                vec3 normal;
                normalGrid(normal, xyz, strides, dim, sds);
                float found_norm = vec3_len(normal);
                // printf("%g\n", found_norm);
                float diff = found_norm - 1.f;
                total += diff * diff;
            }
        }
    }
    long total_elements = dim[0] * dim[1] * dim[2];
    return total / static_cast<float>(total_elements);
}

inline void grid_consistency_loss_diff(
    float *diff_out,
    const long strides[3],
    const long dim[3],
    const float *sds
) {
    extern void __enzyme_autodiff_grid_consistency(
        void *,
        int, const long *,
        int, const long *,
        int, const float *, float *
    );

    __enzyme_autodiff_grid_consistency(
        (void *)grid_consistency_loss,
        enzyme_const, strides,
        enzyme_const, dim,
        enzyme_dup, sds, diff_out
    );
}

inline float sdfPlane(const vec3 pos, const vec3 normal, float height) {
    float dist = vec3_mul_inner(pos, normal) + height;
    return dist;
}

inline float sdfVerticalCapsule(const vec3 origin, float radius, float height) {
    vec3 pos;
    vec3_dup(pos, origin);
    pos[1] -= clamp(pos[1], 0.0f, height);
    return vec3_len(pos) - radius;
}

inline float sdfCylinder(const vec3 pos,float radius, float height) {
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

inline float sdfSphere(const vec3 pos, float radius) {
    return vec3_len(pos) - radius;
}

struct SceneContext {
};

SceneContext *make_scene_context() {
    SceneContext *out = new SceneContext;
    assert(out);
    return out;
}

void free_scene_context(SceneContext *ctx) {
    delete ctx;
}

const long basic_grid_dim = 10;
const long grid_dim[3] = {basic_grid_dim, basic_grid_dim, basic_grid_dim};
const long grid_strides[3] = {basic_grid_dim * basic_grid_dim, basic_grid_dim, 1};

struct SceneParams {
    float object_3_pos[3];
    float color[3];
    float h1, h2;
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
    // float grid[basic_grid_dim][basic_grid_dim][basic_grid_dim];
};

int number_of_scene_params = (int)(sizeof(SceneParams) / sizeof(float));

SceneParams *uninit_scene_params() {
    float nan = nanf("");
    SceneParams *out = new SceneParams;
    for (int p = 0; p < number_of_scene_params; p++) {
        scene_params_set(out, p, nan);
    }
    assert(out);
    return out;
}

void free_scene_params(SceneParams *params) {
    delete params;
}

inline SceneParams *params_from_float_pointer(const float *params) {
    return (SceneParams *)params;
}

const float *float_pointer_from_params(const SceneParams *out) {
    return (float *)out;
}

float *float_pointer_from_params(SceneParams *out) {
    return (float *)out;
}

void scene_params_init(SceneParams *params, const SceneContext *ctx) {
    (void)ctx;
    float *raw_params = float_pointer_from_params(params);
    for (int i = 0; i < number_of_scene_params; i++) {
        raw_params[i] = 0.f;
    }
}

float scene_consistency_loss(const SceneParams *params) {
    (void)params;
    return 0.0;
    // return grid_consistency_loss(grid_strides, grid_dim, &params->grid[0][0][0]);
}

void scene_consistency_gradient(const SceneParams *params, SceneParams *gradient_out) {
    (void)params;
    scene_params_fill(gradient_out, 0.0);
    // grid_consistency_loss_diff(&gradient_out->grid[0][0][0], grid_strides, grid_dim, &params->grid[0][0][0]);
}

void scene_params_copy(SceneParams *out, const SceneParams *params) {
    float *raw_out = float_pointer_from_params(out);
    const float *raw_params = float_pointer_from_params(params);
    for (int i = 0; i < number_of_scene_params; i++) {
        raw_out[i] = raw_params[i];
    }
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

SceneParamsPerChannel *make_ppc(void) {
    SceneParamsPerChannel *ppc = new SceneParamsPerChannel;
    for (int c = 0; c < 3; c++) {
        ppc->rgb[c] = uninit_scene_params();
    }
    return ppc;
}

void set_ppc(SceneParamsPerChannel *ppc, long param, vec3 value) {
    for (int ch = 0; ch < 3; ch++) {
        scene_params_set(ppc->rgb[ch], param, value[ch]);
    }
}

void free_ppc(SceneParamsPerChannel *ppc) {
    for (int ch = 0; ch < 3; ch++) {
        free_scene_params(ppc->rgb[ch]);
    }
    delete ppc;
}

typedef struct {
    float ambient[3];
    float diffuse[3];
    float specular[3];
    float emissive[3];
    bool isReflected;
    float shininess;
} SdfResult;

inline void default_scene_sample(SdfResult *s) {
    vec3_set(s->ambient, 0.0);
    vec3_set(s->diffuse, 0.0);
    vec3_set(s->specular, 0.0);
    vec3_set(s->emissive, 0.0);
    s->isReflected = false;
    s->shininess = 1.0;
}


inline float object_cylinder_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)params;
    (void)ctx;
    vec3 world = {0.6f, 0.2f, -0.2f};
    vec3 local;
    vec3_sub(local, pos, world);
    return sdfCylinder(local, 0.3f, 0.8f);
}
inline void object_cylinder_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    default_scene_sample(sample);
    vec3_set(sample->diffuse, 0.4860f, 0.6310f, 0.6630f);
    vec3_set(sample->ambient, 0.4860f, 0.6310f, 0.6630f);
    vec3_set(sample->specular, 0.8f, 0.8f, 0.8f);
}

inline float object_capsule_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)ctx;
    vec3 world = {-0.4f, 0.3f, 0.5f};
    vec3_add(world, world, params->object_3_pos);
    vec3 local;
    vec3_sub(local, pos, world);
    return sdfVerticalCapsule(local, 0.3f, 0.5f);
}
inline void object_capsule_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    default_scene_sample(sample);
    vec3_set(sample->diffuse, 0.4860f, 0.6310f, 0.6630f);
    vec3_set(sample->ambient, 0.4860f, 0.6310f, 0.6630f);
    vec3_set(sample->specular, 0.8f, 0.8f, 0.8f);
}

// Back:
inline float object_backwall_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)params;
    (void)ctx;
    vec3 plane_normal = {0.0, 0.0, 1.0};
    return sdfPlane(pos, plane_normal, 1.0f);
}
inline void object_backwall_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    default_scene_sample(sample);
    vec3_set(sample->ambient, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->diffuse, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->specular, 0.4f);
}

// Ceiling:
inline float object_topwall_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)params;
    (void)ctx;
    vec3 plane_normal = {0.0, 1.0, 0.0};
    return sdfPlane(pos, plane_normal, 1.0f);
}
inline void object_topwall_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    default_scene_sample(sample);
    vec3_set(sample->ambient, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->diffuse, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->specular, 0.4f);
}

// Left:
inline float object_leftwall_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)params;
    (void)ctx;
    vec3 plane_normal = {1.0, 0.0, 0.0};
    return sdfPlane(pos, plane_normal, 1.0f);
}
inline void object_leftwall_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    default_scene_sample(sample);
    vec3_set(sample->ambient, 0.63f, 0.065f, 0.05f);
    vec3_set(sample->diffuse, 0.63f, 0.065f, 0.05f);
    vec3_set(sample->specular, 0.4f);
}

// Right:
inline float object_rightwall_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)params;
    (void)ctx;
    vec3 plane_normal = {-1.0, 0.0, 0.0};
    return sdfPlane(pos, plane_normal, 1.0f);
}
inline void object_rightwall_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    default_scene_sample(sample);
    vec3_set(sample->ambient, 0.14f, 0.45f, 0.091f);
    vec3_set(sample->diffuse, 0.14f, 0.45f, 0.091f);
    vec3_set(sample->specular, 0.4f);
}

// Floor:
inline float object_bottomwall_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)pos;
    (void)params;
    (void)ctx;
    vec3 plane_normal = {0.0, -1.0, 0.0};
    return sdfPlane(pos, plane_normal, 1.0f);
}
inline void object_bottomwall_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    default_scene_sample(sample);
    vec3_set(sample->ambient, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->diffuse, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->specular, 0.4f);
}

// Grid
inline float object_grid_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)ctx;
    (void)params;
    (void)pos;
    // return sdfGrid(pos, grid_strides, grid_dim, &params->grid[0][0][0]);
    return INFINITY;
}

inline void object_grid_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    default_scene_sample(sample);
    vec3_set(sample->ambient, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->diffuse, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->specular, 0.4f);
}

extern float __enzyme_autodiff_normal(void *, int, const float *, float *, int, const SceneParams *, int, const SceneContext *ctx);

#define SCENE\
    OBJ(object_cylinder_sd,  object_cylinder_mat,  __enzyme_autodiff_normal(\
        (void*)object_cylinder_sd,  enzyme_dup, pos, normal,\
        enzyme_const, params, enzyme_const, ctx))\
    OBJ(object_capsule_sd,   object_capsule_mat,   __enzyme_autodiff_normal(\
        (void*)object_capsule_sd,   enzyme_dup, pos, normal,\
        enzyme_const, params, enzyme_const, ctx))\
    OBJ(object_backwall_sd,  object_backwall_mat,  __enzyme_autodiff_normal(\
        (void*)object_backwall_sd,  enzyme_dup, pos, normal,\
        enzyme_const, params, enzyme_const, ctx))\
    OBJ(object_topwall_sd,   object_topwall_mat,   __enzyme_autodiff_normal(\
        (void*)object_topwall_sd,   enzyme_dup, pos, normal,\
        enzyme_const, params, enzyme_const, ctx))\
    OBJ(object_leftwall_sd,  object_leftwall_mat,  __enzyme_autodiff_normal(\
        (void*)object_leftwall_sd,  enzyme_dup, pos, normal,\
        enzyme_const, params, enzyme_const, ctx))\
    OBJ(object_rightwall_sd, object_rightwall_mat, __enzyme_autodiff_normal(\
        (void*)object_rightwall_sd, enzyme_dup, pos, normal,\
        enzyme_const, params, enzyme_const, ctx))\
    OBJ(object_bottomwall_sd,object_bottomwall_mat,__enzyme_autodiff_normal(\
        (void*)object_bottomwall_sd,enzyme_dup, pos, normal,\
        enzyme_const, params, enzyme_const, ctx))
    // OBJ(object_grid_sd,object_grid_mat,__enzyme_autodiff_normal(\
    //     (void*)object_grid_sd,enzyme_dup, pos, normal,\
    //     enzyme_const, params, enzyme_const, ctx))

inline float scene_sample_sdf(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    float ret = INFINITY;
    #define OBJ(sd, mat, norm) do { ret = fminf(ret, sd(pos, params, ctx)); } while (0);
    SCENE
    #undef OBJ
    return ret;
}

inline float get_normal_from(vec3 out_normal, const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    float ret = INFINITY;
    #define OBJ(sd, mat, norm) do {\
        vec3 normal; \
        vec3_set(normal, 0.0); \
        float distance = sd(pos, params, ctx);\
        norm;\
        if (distance < ret) { \
            vec3_dup(out_normal, normal); \
            ret = distance;\
        } \
    } while (0);
    SCENE
    #undef OBJ
    return ret;
}

/** writes the output in the first paramemter */
inline float compose_scene_sample(SdfResult *a, const SdfResult *b, float distance_a, float distance_b) {
    if (distance_a < distance_b) {
        return distance_a; // nothing to be done
    }
    vec3_dup(a->ambient, b->ambient);
    vec3_dup(a->diffuse, b->diffuse);
    vec3_dup(a->specular, b->specular);
    a->shininess = b->shininess;
    return distance_b;
}

inline void scene_sample(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *res) {
    SdfResult working;
    default_scene_sample(&working);
    float res_distance = INFINITY;
    #define OBJ(sd, mat, norm) do { \
        float distance_a = sd(pos, params, ctx);\
        mat(pos, params, ctx, &working); \
        res_distance = compose_scene_sample(res, &working, res_distance, distance_a); \
    } while(0);
    SCENE
    #undef OBJ
}

inline float directional_derivative(const vec3 origin, const vec3 direction, float t, const SceneParams *params, const SceneContext *ctx) {
    vec3 sample_point;
    ray_step(sample_point, origin, direction, t);
    vec3 normal;
    get_normal_from(normal, sample_point, params, ctx);
    return vec3_mul_inner(direction, normal);
}

typedef struct {
    vec3 origin, direction;
    const SceneParams *params;
    const SceneContext *ctx;
    int iter_cap;
} CriticalPointContext;

inline float critical_point_bisect_midpoint(float a, float b, void *context) {
    CriticalPointContext *ctx = (CriticalPointContext *)context;
    float dir_a = directional_derivative(ctx->origin, ctx->direction, a, ctx->params, ctx->ctx);
    float dir_b = directional_derivative(ctx->origin, ctx->direction, b, ctx->params, ctx->ctx);
    // have points (a, dir_a) and (b, dir_b). want to find the zero crossing.
    float denom = dir_a - dir_b;
    if (fabsf(denom) < 1e-5) {
        return a; // arbitrarily choose one of the endpoints
    }
    float numer = dir_a * b - a * dir_b;
    return numer / denom;
}

inline BisectAffinity critical_point_bisect_affinity(float t, int iter_count, void *context) {
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
        if (scene_sample_sdf(pos, ctx->params, ctx->ctx) > preliminary_distance_threshold) {
            return BISECT_STOP;
        }
    }
    float dir_t = directional_derivative(ctx->origin, ctx->direction, t, ctx->params, ctx->ctx);
    bool is_approaching = dir_t < 0;
    return is_approaching ? BISECT_RIGHT : BISECT_LEFT;
}

typedef struct {
    bool found_critical_point;
    float t_if_found_critical_point;
} SearchResult;

void search_for_critical_point(
    const vec3 origin, const vec3 direction, const SceneParams *params, const SceneContext *ctx, float t_min, float t_max, SearchResult *ret
) {
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
    if (scene_sample_sdf(best_pos, params, ctx) > distance_threshold) {
        return;
    }
    float best_dir_t = directional_derivative(origin, direction, best_t, params, ctx);
    if (fabsf(best_dir_t) > directional_derivative_threshold) {
        return;
    }
    ret->found_critical_point = true;
    ret->t_if_found_critical_point = best_t;
}

float sdf_theta_wrapper(const vec3 pos, const float *params, const SceneContext *ctx) {
    SceneParams *scene_params = params_from_float_pointer(params);
    return scene_sample_sdf(pos, scene_params, ctx);
}

extern void __enzyme_autodiff_theta(void *, int,const float *, int, const float *, float *, int, const SceneContext *ctx);

void diff_sdf(const vec3 pos, SceneParams *paramsOut, const SceneParams *paramsIn, const SceneContext *ctx) {
    const float *raw_params = float_pointer_from_params(paramsIn);
    float* draw_params = float_pointer_from_params(paramsOut);
    __enzyme_autodiff_theta(
        (void*)sdf_theta_wrapper,
        enzyme_const, pos,
        enzyme_dup, raw_params, draw_params,
        enzyme_const, ctx
    );
}

inline void phongLight(vec3 radiance, const vec3 looking, const vec3 normal, const SdfResult *sample, const SceneParams *params) {
    float lightColors[3][3] = {
        {.2f, .2f, .2f},
        {.2f, .2f, .2f},
        {.2f, .2f, .2f},
    };
    vec3_add(lightColors[0], lightColors[0], params->color);
    // world[0] -= params->object_3_pos[0];
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

inline IntersectionResult trace_ray_get_intersection(
    const vec3 origin,
    const vec3 direction,
    const SceneParams *params,
    const SceneContext *ctx
) {
    float t = 0.0;
    IntersectionResult ret = { false, 0.0 };
    for (int i = 0; i < number_of_steps; i++) {
        vec3 pos;
        ray_step(pos, origin, direction, t);
        float distance = scene_sample_sdf(pos, params, ctx);
        if (distance < contact_threshold) {
            ret.found_intersection = true;
            ret.intersection_t = t;
            break;
        }
        t += distance;
    }
    return ret;
}

IntersectionResult trace_ray_get_critical_point(
    SearchResult *critical_point,
    const vec3 origin,
    const vec3 direction,
    const SceneParams *params,
    const SceneContext *ctx
) {
    float t = 0.0;
    float previous_t = 0.0;
    critical_point->found_critical_point = false;
    IntersectionResult ret = { false, 0.0 };
    for (int i = 0; i < number_of_steps; i++) {
        float dir_previous_t = directional_derivative(origin, direction, previous_t, params, ctx);
        float dir_t = directional_derivative(origin, direction, t, params, ctx);
        // look for a sign change in the directional derivative
        if (!critical_point->found_critical_point && ((dir_previous_t < 0) && (dir_t > 0))) {
            // let's try to find critical_point between t and previous_t;
            search_for_critical_point(origin, direction, params, ctx, previous_t, t, critical_point);
        }
        vec3 pos;
        ray_step(pos, origin, direction, t);
        float distance = scene_sample_sdf(pos, params, ctx);
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

void sample_points_on_plane(int num_points, const vec3 center_on_plane, float stddev, vec3 *output_points, RandomState* random) {
    if (fabsf(center_on_plane[1] - 1.0f) > 1e-6) {
         fprintf(stderr, "Warning: Center point provided is not on the plane y=1.0\n");
         // Adjust center y to be exactly on the plane if desired, or proceed.
         // center_on_plane.y = 1.0f; // Optional correction
    }

    for (int i = 0; i < num_points; ++i) {
        // Generate standard normal random numbers (mean 0, stddev 1)
        float rand_x_std = generate_normal_random(random);
        float rand_z_std = generate_normal_random(random);

        // Scale by standard deviation and add to the center coordinates
        float point_x = center_on_plane[0] + stddev * rand_x_std;
        float point_z = center_on_plane[2] + stddev * rand_z_std;
        float point_y = center_on_plane[1]; // Keep y fixed on the plane

        vec3_set(output_points[i], point_x, point_y, point_z);
    }
}

float normal_pdf(float x, float mu, float sigma) {
    if (sigma <= 0.0f) { return 0.0f; } // Avoid division by zero / invalid input
    float diff = x - mu;
    float exponent = -(diff * diff) / (2.0f * sigma * sigma);
    float scaling_factor = 1.0f / (sigma * sqrtf(2.0f * lm_pi));
    return scaling_factor * expf(exponent);
}
// --------------------------------------------------------------------

/**
 * @brief Calculates the joint PDF for a point on a plane, assuming its X and Z
 * coordinates were sampled independently from N(center.x, sigma^2)
 * and N(center.z, sigma^2) respectively.
 * Assumes point.y is fixed and matches center.y implicitly.
 *
 * @param point The vec3 point (e.g., {px, py, pz}) to evaluate.
 * @param center The vec3 center of the distribution on the plane (e.g., {cx, cy, cz}).
 * @param sigma The standard deviation used for sampling X and Z.
 * @return float The joint probability density value for the XZ coordinates.
 */
float plane_pdf_independent_xz(const vec3 point, const vec3 center, float sigma) {
    // Assumes point and center are float[3] or compatible struct access
    // Use point[0] or point.x depending on your vec3 definition
    float pdf_x = normal_pdf(point[0], center[0], sigma);
    float pdf_z = normal_pdf(point[2], center[2], sigma);

    // Joint PDF is the product due to independence
    return pdf_x * pdf_z;
}

void computeDirectLighting(vec3 hitPoint, vec3 hitNormal, SdfResult sample,const SceneParams *params, const SceneContext *ctx, vec3 directLighting, RandomState* random){
    const int num_samples = 1;
    vec3 sampled_points[num_samples];

    vec3 center;
    vec3_set(center, 0.0f, 1.0f, 0.0f); // Center of distribution on the plane y=1.0
    float standard_deviation = 0.5f;     // Adjust spread as needed

    sample_points_on_plane(num_samples, center, standard_deviation, sampled_points,random);

    vec3 contribution;
    vec3_set(contribution,0.0f);

    for(int i = 0; i < num_samples;i++){
        vec3 currentSamplePoint;

        vec3_dup(currentSamplePoint,sampled_points[i]);

        vec3 light_direction;

        vec3_sub(light_direction,currentSamplePoint,hitPoint);

        float distance2 = vec3_mul_inner(light_direction, light_direction);

        vec3_norm(light_direction,light_direction);

        IntersectionResult shadow_hit_point = trace_ray_get_intersection(hitPoint,light_direction,params, ctx);
        if(shadow_hit_point.found_intersection){
            const float ray_epsilon1 = 1e-4f;

            vec3 shadow_hit_intersection_point;
            ray_step(shadow_hit_intersection_point, hitPoint, light_direction, shadow_hit_point.intersection_t);

            vec3 shawdow_normal;
            get_normal_from(shawdow_normal,shadow_hit_intersection_point,params, ctx);

            vec3 offset_position;
            vec3_scale(offset_position, shawdow_normal, ray_epsilon1);
            vec3_add(shadow_hit_intersection_point,shadow_hit_intersection_point,offset_position);

            vec3 diff_hit;
            vec3_sub(diff_hit,shadow_hit_intersection_point,hitPoint);
            float hitDistance = vec3_len(diff_hit);

            vec3 diff_light;
            vec3_sub(diff_light,currentSamplePoint,hitPoint);
            float distToLight = vec3_len(diff_light);

            if(fabsf(hitDistance - distToLight) < 1e-4f){

                float costheta = fmaxf(0.f, vec3_mul_inner(light_direction,hitNormal));

                vec3 emission;
                vec3_set(emission, 0.9f);
                vec3 lightNormal;
                vec3_set(lightNormal,0.0f,1.0f,0.0f);

                float cosThetaPrime = fmax(0.0f, vec3_mul_inner(light_direction, lightNormal));

                vec3 bsdf;
                vec3_set(bsdf,1.0f);

                vec3_scale(bsdf,sample.diffuse, 1.0f / lm_pi);

                float pdf = plane_pdf_independent_xz(currentSamplePoint,center,standard_deviation);

                float last_term = costheta * cosThetaPrime / ((distance2 + 1e-6f) * pdf);

                vec3_scale(bsdf, bsdf,last_term);
                vec3_cwiseProduct(bsdf,bsdf,emission);

                vec3_add(contribution,contribution,bsdf);
            }
        }
    }
    vec3_scale(contribution,contribution,1.0f/num_samples);
    vec3_dup(directLighting,contribution);
}

// void toneMap(vec3 tone_map_color, const vec3 color) {
//     vec3 map_weights = {0.2126f, 0.7152f, 0.0722f};
//     float luminance = vec3_mul_inner(color, map_weights);
//     // Reinhard Tone Mapping
//     float toneMappedLuminance = luminance / (1.0f + luminance);
//     // Scale the color channels by the tone-mapped luminance
//     if (luminance > 1e-6f) { // Avoid division by zero
//         float scale = toneMappedLuminance / luminance;
//         vec3_scale(tone_map_color, color, scale);
//     } else {
//         vec3_dup(tone_map_color, color);
//     }
// }

void gammaCorrect(vec3 gamma_correct_color, const vec3 color) {
    float gamma = 2.2f;
    gamma_correct_color[0] = powf(fmaxf(0.0f, color[0]), 1.0f / gamma);
    gamma_correct_color[1] = powf(fmaxf(0.0f, color[1]), 1.0f / gamma);
    gamma_correct_color[2] = powf(fmaxf(0.0f, color[2]), 1.0f / gamma);
}

inline void get_radiance_phong(
    vec3 radiance,
    const IntersectionResult *intersection,
    const vec3 origin,
    const vec3 direction,
    const SceneParams *params,
    const SceneContext *ctx
) {
    if (!intersection->found_intersection) {
        vec3_set(radiance, 0.f);
        return;
    }
    vec3 current_position;
    ray_step(current_position, origin, direction, intersection->intersection_t);
    SdfResult sample;
    scene_sample(current_position, params, ctx, &sample);
    vec3 normal;
    get_normal_from(normal, current_position, params, ctx);
    phongLight(radiance, direction, normal, &sample, params);
}

void params_nan_to_num(SceneParamsPerChannel *ppc, float num) {
    for (int ch = 0; ch < 3; ch++) {
        for (long p = 0; p < number_of_scene_params; p++) {
            float param = scene_parameter_get(ppc->rgb[ch], p);
            if (isnan(param)) {
                scene_params_set(ppc->rgb[ch], p, num);
            }
        }
    }
}

GradientImage make_gradient_image(long image_width, long image_height) {
    const long num_subpixels = 3;
    long num_floats = number_of_scene_params * image_width * image_height * num_subpixels;
    assert(num_floats > 0);
    float *buf = new float[(size_t)num_floats];
    assert(buf);
    for (long i = 0; i < num_floats; i++) {
        buf[i] = 0.f;
    }
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
    SceneParamsPerChannel *ppc = make_ppc();
    for (long r = 0; r < image->image_height; r++) {
        for (long c = 0; c < image->image_width; c++) {
            gradient_image_get(ppc, gradient, r, c);
            vec3 radiance;
            for (int ch = 0; ch < 3; ch++) {
                float param_value = scene_parameter_get(ppc->rgb[ch], parameter_no);
                radiance[ch] = param_value;
            }
            vec3 half;
            vec3_set(half, 0.5f);
            vec3_scale(radiance, radiance, 0.5f);
            vec3_add(radiance, radiance, half);
            image_set(image, r, c, radiance);
        }
    }
    free_ppc(ppc);
}

void free_gradient_image(GradientImage *image) {
    delete[] image->buf;
}

inline long gradient_image_get_index(const GradientStrides *s, long r, long c, long subpixel, long param) {
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

float render_get_radiance_wrapper(
    const IntersectionResult *intersection,
    const vec3 origin,
    const vec3 direction,
    const float *raw_params,
    const SceneContext *ctx,
    int ch
) {
    const SceneParams *params = params_from_float_pointer(raw_params);
    vec3 naive_intersection;
    ray_step(naive_intersection, origin, direction, intersection->intersection_t);
    float sdf_at_intersection = scene_sample_sdf(naive_intersection, params, ctx);
    float dsdfdt = directional_derivative(origin, direction, intersection->intersection_t, params, ctx);
    float first_order_update = - sdf_at_intersection / dsdfdt;
    IntersectionResult sensitivity_adjusted;
    sensitivity_adjusted.found_intersection = intersection->found_intersection;
    sensitivity_adjusted.intersection_t = intersection->intersection_t + first_order_update;
    vec3 radiance;
    get_radiance_phong(radiance, &sensitivity_adjusted, origin, direction, params, ctx);
    return radiance[ch];
}

void render_pixel_phong(
    vec3 real,
    const vec3 origin,
    const vec3 direction,
    const SceneParams *params,
    const SceneContext *ctx
) {
    IntersectionResult intersection = trace_ray_get_intersection(origin, direction, params, ctx);
    get_radiance_phong(real, &intersection, origin, direction, params, ctx);
}

extern void __enzyme_autodiff_radiance(
    void *,
    int, const IntersectionResult *,
    int, const float *,
    int, const float *,
    int, const float *, float *,
    int, const SceneContext *ctx,
    int, int
);

typedef struct {
    vec3 pos;
    vec3 dir;
    vec3 contribution;
} Segment;

void step_segment(Segment &segment, float t) {
    ray_step(segment.pos, segment.pos, segment.dir, t);
}

std::vector<Segment> getSecondaryPath(
    const vec3 origin,
    const vec3 direction,
    RandomState *random,
    const SceneParams *params,
    const SceneContext *ctx
) {
    std::vector<Segment> path;

    vec3 current_position;
    vec3_dup(current_position, origin);
    vec3 current_direction;
    vec3_dup(current_direction,direction);
    vec3 previous_direction;
    vec3_set(previous_direction,0.0f);

    int iter = 0;

    while(true){
        iter++;
        // printf("%d\n", iter);
        IntersectionResult hitPoint = trace_ray_get_intersection(current_position,current_direction,params, ctx);
        if(!hitPoint.found_intersection){
            return path;
        }
        //store current ray
        Segment newSegment;
        vec3_dup(newSegment.pos,current_position);
        vec3_dup(newSegment.dir,current_direction);

        vec3 directLighting;

        if(iter == 0){
            vec3_set(directLighting, 0.0f);
        }else{
            SdfResult hit_sample;
            scene_sample(current_position, params, ctx, &hit_sample);

            vec3 hit_normal;
            get_normal_from(hit_normal, current_position, params, ctx);

            computeDirectLighting(current_position,hit_normal,hit_sample,params,ctx, directLighting,random);

        }
        vec3_dup(newSegment.contribution,directLighting);

        path.push_back(newSegment);

        if(random_next_float(random)> pathContinuationProb){
            break;
        }

        //find intersection point

        vec3 hit_intersection_point;
        ray_step(hit_intersection_point, current_position, current_direction, hitPoint.intersection_t);

        SdfResult sample;
        scene_sample(hit_intersection_point, params, ctx, &sample);

        vec3 normal;
        get_normal_from(normal, hit_intersection_point, params, ctx);


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
        vec3_dup(previous_direction,current_direction);
        vec3_dup(current_direction,newDir);
    }

    return path;
}

void get_radiance_tracing(
    vec3 radiance,
    const std::vector<Segment> &path,
    const SceneParams *params,
    const SceneContext *ctx
) {
    vec3 spectralFilter;
    vec3_set(spectralFilter, 1.f);
    vec3 intensity;
    vec3_set(intensity, 0.f);
    for(size_t i = 0; i < path.size(); i++){
        if(i != path.size() - 1){
            vec3 hitPosition;
            vec3_dup(hitPosition, path[i+1].pos);
            vec3 hitDirection;
            vec3_dup(hitDirection,path[i].dir);
            vec3 wi;
            vec3_dup(wi,path[i+1].dir);

            SdfResult sample;
            scene_sample(hitPosition, params, ctx, &sample);
            vec3 normal;
            get_normal_from(normal, hitPosition, params, ctx);

            vec3 emissive;
            vec3_dup(emissive, sample.emissive);

            vec3 emissive_part;
            vec3_cwiseProduct(emissive_part, emissive, spectralFilter);

            vec3_add(intensity, intensity, emissive_part);

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
}

float render_pixel_tracing_wrapper(
    const std::vector<Segment> &path,
    const float *raw_params,
    const SceneContext *ctx,
    int ch
) {
    const SceneParams *params = params_from_float_pointer(raw_params);
    vec3 radiance;
    get_radiance_tracing(radiance, path, params, ctx);
    return radiance[ch];
}

void render_gradient_phong(
    vec3 real,
    const vec3 origin,
    const vec3 direction,
    SceneParamsPerChannel *params_per_channel,
    const SceneParams *params,
    const SceneContext *ctx,
    RandomState *rng
) {
    (void)rng;
    SearchResult critical_point;
    IntersectionResult intersection = trace_ray_get_critical_point(&critical_point, origin, direction, params, ctx);
    const float *raw_params = float_pointer_from_params(params);
    for (int ch = 0; ch < 3; ch++) {
        SceneParams *fill = params_per_channel->rgb[ch];
        scene_params_fill(fill, 0.f);
        float *raw_fill = float_pointer_from_params(fill);
        extern void __enzyme_autodiff_radiance(
            void *,
            int, const IntersectionResult *,
            int, const float *,
            int, const float *,
            int, const float *, float *,
            int, const SceneContext *ctx,
            int, int
        );
        __enzyme_autodiff_radiance(
            (void*)render_get_radiance_wrapper,
            enzyme_const, &intersection,
            enzyme_const, origin,
            enzyme_const, direction,
            enzyme_dupnoneed, raw_params, raw_fill,
            enzyme_const, ctx,
            enzyme_const, ch
        );
    }
    get_radiance_phong(real, &intersection, origin, direction, params, ctx);
    if (critical_point.found_critical_point) {
        SceneParams *dummy_params = uninit_scene_params();
        scene_params_fill(dummy_params, 0.f);
        vec3 y_star;
        ray_step(y_star, origin, direction, critical_point.t_if_found_critical_point);
        SdfResult sample;
        scene_sample(y_star, params, ctx, &sample);
        vec3 normal;
        get_normal_from(normal, y_star, params, ctx);
        vec3 y_star_radiance;
        phongLight(y_star_radiance, direction, normal, &sample, params);
        diff_sdf(y_star, dummy_params, params, ctx);
        vec3 deltaL;
        vec3_sub(deltaL, y_star_radiance, real);
        vec3_scale(deltaL, deltaL, -1.f / distance_threshold);
        outer_product_add_assign(params_per_channel, dummy_params, deltaL);
        free_scene_params(dummy_params);
    }
    params_nan_to_num(params_per_channel, 0.0);
}

void render_pixel_tracing(
    vec3 real,
    const vec3 origin,
    const vec3 direction,
    SceneParamsPerChannel *params_per_channel,
    const SceneParams *params,
    const SceneContext *ctx,
    RandomState *rng
) {
    SearchResult critical_point;
    IntersectionResult intersection = trace_ray_get_critical_point(&critical_point, origin, direction, params, ctx);
    (void)intersection;
    std::vector<Segment> path = getSecondaryPath(origin, direction, rng, params, ctx);
    get_radiance_tracing(real, path, params, ctx);
    const float *raw_params = float_pointer_from_params(params);
    for (int ch = 0; ch < 3; ch++) {
        SceneParams *fill = params_per_channel->rgb[ch];
        scene_params_fill(fill, 0.f);
        float *raw_fill = float_pointer_from_params(fill);
        extern void __enzyme_autodiff_radiance2(
            void *,
            int, const std::vector<Segment> &,
            int, const float *, float *,
            int, const SceneContext *ctx,
            int, int
        );
        __enzyme_autodiff_radiance2(
            (void*)render_pixel_tracing_wrapper,
            enzyme_const, path,
            enzyme_dupnoneed, raw_params, raw_fill,
            enzyme_const, ctx,
            enzyme_const, ch
        );
    }
    if (critical_point.found_critical_point && path.size() > 0) {
        SceneParams *dummy_params = uninit_scene_params();
        scene_params_fill(dummy_params, 0.f);
        vec3 y_star;
        ray_step(y_star, origin, direction, critical_point.t_if_found_critical_point);
        vec3 y_star_radiance;
        std::vector<Segment> other_path(path);
        step_segment(other_path.at(0), critical_point.t_if_found_critical_point);
        get_radiance_tracing(y_star_radiance, other_path, params, ctx);
        diff_sdf(y_star, dummy_params, params, ctx);
        vec3 deltaL;
        vec3_sub(deltaL, y_star_radiance, real);
        vec3_scale(deltaL, deltaL, -1.f / distance_threshold);
        outer_product_add_assign(params_per_channel, dummy_params, deltaL);
        free_scene_params(dummy_params);
    }
    params_nan_to_num(params_per_channel, 0.0);
}

inline long get_index(const Strides *s, long r, long c, long p) {
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
    long found_width;
    if (fscanf(f, "%ld", &found_width) != 1) return 0;
    assert(found_width == image->image_width);

    skip_whitespace_and_comments(f);
    long found_height;
    if (fscanf(f, "%ld", &found_height) != 1) return 0;
    assert(found_height == image->image_height);

    skip_whitespace_and_comments(f);
    int maxval;
    if (fscanf(f, "%d", &maxval) != 1 || maxval != 255) {
        fprintf(stderr, "Unsupported max color value (only 255 supported)\n");
        return 0;
    }
    // Skip the single whitespace character after maxval
    fgetc(f);
    uint8_t *buf = new uint8_t[(size_t)image->num_floats];
    size_t bytes_read = fread(buf, sizeof(uint8_t), (size_t)image->num_floats, f);
    assert(bytes_read == (size_t)image->num_floats);
    for (long i = 0; i < image->num_floats; i++) {
        image->buf[i] = (float)buf[i] * (1.f / 255.f);
    }
    delete[] buf;
    return 1; // success
}
// End PPM Parser from ChatGPT ///////
//////////////////////////////////////

Image make_image(long image_width, long image_height) {
    Image image;
    long num_floats = image_width * image_height * 3;
    assert(num_floats > 0);
    float *buf = new float[(size_t)num_floats];
    assert(buf);
    Strides strides;
    strides.row_stride = image_width * 3;
    strides.col_stride = 3;
    strides.subpixel_stride = 1;
    image.strides = strides;
    image.image_width = image_width;
    image.image_height = image_height;
    image.buf = buf;
    image.num_floats = num_floats;
    return image;
}

void free_image(Image *image) {
    delete[] image->buf;
}

// https://nullprogram.com/blog/2017/11/03/
void image_write_ppm(const Image *image, FILE *f) {
    assert(f);
    fprintf(f, "P6\n%ld %ld\n255\n", image->image_width, image->image_height);
    assert(image->num_floats > 0);
    uint8_t *buf = new uint8_t[(size_t)image->num_floats];
    for (long i = 0; i < image->num_floats; i++) {
        uint8_t value = (uint8_t)(clamp(roundf(image->buf[i] * 255.f), 0.f, 255.f));
        buf[i] = value;
    }
    fwrite(buf, 1, (size_t)image->num_floats, f);
    delete[] buf;
    fflush(f);
}

void image_set(Image *image, long ir, long ic, const vec3 radiance) {
    assert(ir >= 0);
    assert(ir < image->image_height);
    assert(ic >= 0);
    assert(ic < image->image_width);
    for (long p = 0; p < 3; p++) {
        long index = get_index(&image->strides, ir, ic, p);
        image->buf[index] = radiance[p];
    }
}

void image_get(vec3 radiance, const Image *image, long ir, long ic) {
    assert(ir >= 0);
    assert(ir < image->image_height);
    assert(ic >= 0);
    assert(ic < image->image_width);
    for (long p = 0; p < 3; p++) {
        long index = get_index(&image->strides, ir, ic, p);
        radiance[p] = image->buf[index];
    }
}

void image_copy(const Image *src, Image *dst) {
    for (long ir = 0; ir < src->image_height; ir++) {
        for (long ic = 0; ic < src->image_width; ic++) {
            vec3 radiance;
            image_get(radiance, src, ir, ic);
            image_set(dst, ir, ic, radiance);
        }
    }
}

struct Renderer {
    float world_from_camera[4][4];
    vec3 camera_position;
    long image_width, image_height;
};

Renderer *make_renderer(long image_width, long image_height) {
    float aspect = (float)image_width / (float)image_height ;
    float near_clip = 0.1f;
    float far_clip = 100.0f;
    float y_fov = 0.785f;
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
    Renderer *ret = new Renderer;
    mat4x4_dup(ret->world_from_camera, world_from_camera);
    vec3_dup(ret->camera_position, camera_position);
    ret->image_width = image_width;
    ret->image_height = image_height;
    return ret;
}

void free_renderer(Renderer *renderer) {
    delete renderer;
}

struct PixelRenderer {
    float camera_position[3];
    float direction[3];
};

PixelRenderer make_pixel_renderer(const Renderer *renderer, long ir, long ic) {
    float r = (float)ir;
    float c = (float)ic;
    float device_x = lerp(c, 0.0, (float)renderer->image_width, -1.0, 1.0);
    float device_y = lerp(r, 0.0, (float)renderer->image_height, -1.0, 1.0);
    vec4 unprojected = {device_x, device_y, 1.0, 1.0};
    vec4 homo;
    mat4x4_mul_vec4(homo, renderer->world_from_camera, unprojected);
    vec3 far;
    dehomogenize(far, homo);
    vec3 direction;
    vec3_sub(direction, far, renderer->camera_position);
    vec3_norm(direction, direction);
    PixelRenderer ret;
    vec3_dup(ret.camera_position, renderer->camera_position);
    vec3_dup(ret.direction, direction);
    return ret;
}

void render_image_phong(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *rng) {
    (void)rng;
    Renderer *renderer = make_renderer(real->image_width, real->image_height);
    #pragma omp parallel for
    for (long ir = 0; ir < real->image_height; ir++) {
        SceneParamsPerChannel *ppc = make_ppc();
        for (long ic = 0; ic < real->image_width; ic++) {
            PixelRenderer pr = make_pixel_renderer(renderer, ir, ic);
            vec3 out_real;
            render_gradient_phong(out_real, pr.camera_position, pr.direction, ppc, params, ctx, rng);
            image_set(real, ir, ic, out_real);
            gradient_image_set(ppc, gradient, ir, ic);
        }
        free_ppc(ppc);
    }
    free_renderer(renderer);
}

void render_image_tracing(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *rng) {
    (void)rng;
    Renderer *renderer = make_renderer(real->image_width, real->image_height);
    #pragma omp parallel for
    for (long ir = 0; ir < real->image_height; ir++) {
        SceneParamsPerChannel *ppc = make_ppc();
        for (long ic = 0; ic < real->image_width; ic++) {
            PixelRenderer pr = make_pixel_renderer(renderer, ir, ic);
            vec3 out_real;
            render_pixel_tracing(out_real, pr.camera_position, pr.direction, ppc, params, ctx, rng);
            image_set(real, ir, ic, out_real);
            gradient_image_set(ppc, gradient, ir, ic);
        }
        free_ppc(ppc);
    }
    free_renderer(renderer);
}

void free_pixel_renderer(PixelRenderer *renderer) {
    delete renderer;
}

void finite_differences(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *rng) {
    (void)rng;
    Renderer *renderer = make_renderer(real->image_width, real->image_height);

    #pragma omp parallel for
    for (long ir = 0; ir < real->image_height; ir++) {
        SceneParams *working = uninit_scene_params();
        SceneParamsPerChannel *ppc = make_ppc();
        for (long ic = 0; ic < real->image_width; ic++) {
            PixelRenderer pr = make_pixel_renderer(renderer, ir, ic);
            for (long p = 0; p < number_of_scene_params; p++) {
                float param_value = scene_parameter_get(params, p);
                scene_params_copy(working, params);
                scene_params_set(working, p, param_value - finite_difference_epsilon);
                vec3 real_left;
                render_pixel_phong(real_left, pr.camera_position, pr.direction, working, ctx);
                vec3 real_right;
                scene_params_copy(working, params);
                scene_params_set(working, p, param_value + finite_difference_epsilon);
                render_pixel_phong(real_right, pr.camera_position, pr.direction, working, ctx);

                vec3 gradient;
                vec3_sub(gradient, real_right, real_left);
                vec3_scale(gradient, gradient, 1.f / (2.f * finite_difference_epsilon));

                for (int ch = 0; ch < 3; ch++) {
                    scene_params_set(ppc->rgb[ch], p, gradient[ch]);
                }
            }
            gradient_image_set(ppc, gradient, ir, ic);
        }
        free_scene_params(working);
        free_ppc(ppc);
    }
    free_renderer(renderer);
}

void chromaticAberrationFilter(long rShift, long gShift, long bShift, Image* outImage, const Image* real){
    for (long ir = 0; ir < real->image_height; ir++) {
        for (long ic = 0; ic < real->image_width; ic++) {
            vec3 pixel;
            image_get(pixel,real,ir,ic);

            long rShiftedCol = clampl(ic + rShift, 0, real->image_width - 1);
            long gShiftedCol = clampl(ic + gShift, 0, real->image_width - 1);
            long bShiftedCol = clampl(ic + bShift, 0, real->image_width - 1);

            vec3 rChannelPixel;
            vec3 gChannelPixel;
            vec3 bChannelPixel;
            image_get(rChannelPixel, real, ir, rShiftedCol);
            image_get(gChannelPixel, real, ir, gShiftedCol);
            image_get(bChannelPixel, real, ir, bShiftedCol);

            vec3 outputPixel;
            outputPixel[0] = rChannelPixel[0];
            outputPixel[1] = gChannelPixel[1];
            outputPixel[2] = bChannelPixel[2];

            image_set(outImage, ir, ic, outputPixel);
        }
    }
}

float* generateBlurkernel(float sigma){
    float* kernel = nullptr;
    int kernel_size = 2 * blur_radius + 1;
    assert(kernel_size > 0);
    kernel = new float[static_cast<size_t>(kernel_size)];
    float sum = 0.0f;
    // Calculate the kernel values
    for (int i = 0; i < kernel_size; ++i) {
        float x = static_cast<float>(i - 2);  // Adjust x to the correct range
        float exponent = -(x * x) / (2 * sigma * sigma);
        kernel[i] = (1.0f / (sqrt(2 * lm_pi) * sigma)) * exp(exponent);
        sum += kernel[i];  // Accumulate the sum for normalization
    }
    for (int i = 0; i < kernel_size; ++i) {
        kernel[i] /= sum;
    }
    return kernel;
}

void blurFilter(Image* outImage, const Image* real, float sigma){
    float* kernel = generateBlurkernel(sigma);
    long kernel_size = 2 * blur_radius + 1;
    for (long ir = 0; ir < real->image_height; ir++) {
        for (long ic = 0; ic < real->image_width; ic++) {
            float redAcc = 0.0f, greenAcc = 0.0f, blueAcc = 0.0f;
            for (long k = 0; k < kernel_size; ++k) {
                long pixelCol = clampl(ic + k - blur_radius, 0L, static_cast<long>(real->image_width - 1));
                vec3 pixel;
                image_get(pixel,real,ir,pixelCol);
                float blurWeight = kernel[k];
                redAcc += pixel[0] * blurWeight;
                greenAcc += pixel[1]* blurWeight;
                blueAcc += pixel[2] * blurWeight;
            }
            vec3 out;
            out[0] = fmaxf(0.0f, fminf(1.f, redAcc));
            out[1] = fmaxf(0.0f, fminf(1.f, greenAcc));
            out[2] = fmaxf(0.0f, fminf(1.f, blueAcc));
            image_set(outImage, ir, ic, out);
        }
    }
    delete[] kernel;
}

void render_image_effects(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *rng) {
    (void)gradient;
    Image scratch = make_image(real->image_width, real->image_height);
    render_image_phong(real, gradient, params, ctx, rng);
    // chromaticAberrationFilter(5, 2, 10, &scratch, real);
    float sigma = 1.5f;
    blurFilter(&scratch, real, sigma);
    image_copy(&scratch, real);
    free_image(&scratch);
}
