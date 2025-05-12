#include <cassert>
#include <functional>

#include "linmath.h"
#include "util.h"
#include "params.h"
#include "scene.h"
#include "linmath.h"

inline void default_scene_sample(SdfResult *s) {
    vec3_set(s->ambient, 0.0);
    vec3_set(s->diffuse, 0.0);
    vec3_set(s->specular, 0.0);
    vec3_set(s->emissive, 0.0);
    s->isReflected = false;
    s->shininess = 1.0;
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

float sdfSquare(const vec3 pos, const vec3 center) {
    vec3 displaced_pos;
    vec3_sub(displaced_pos,pos,center);
    // Distance to the sides of the square in the xz-plane
    float halfSide = 0.5f;
    float halfThickness = 0.5f / 2.0f;

    // Distance to the sides of the square in the xz-plane
    float distToRight = displaced_pos[0] - halfSide;
    float distToLeft = -displaced_pos[0] - halfSide;
    float distToFront = displaced_pos[2];
    float distToBack = -displaced_pos[2];

    // Distance to the top and bottom planes (along the y-axis)
    float distToTop = displaced_pos[1] - halfThickness;
    float distToBottom = -displaced_pos[1] - halfThickness;

    // The SDF is the maximum of these distances
    return fmaxf(fmaxf(fmaxf(distToRight, distToLeft), fmaxf(distToFront, distToBack)), fmaxf(distToTop, distToBottom));
}

float sdfTriangle(const vec3 pos, const vec3 a, const vec3 b, const vec3 c) {
    vec3 ba;
    vec3_sub(ba, b, a);
    vec3 pa;
    vec3_sub(pa, pos, a);
    vec3 cb;
    vec3_sub(cb, c, b);
    vec3 pb;
    vec3_sub(pb, pos, b);
    vec3 ac;
    vec3_sub(ac, a, c);
    vec3 pc;
    vec3_sub(pc, pos, c);
    vec3 nor;
    vec3_mul_cross(nor, ba, ac);
    vec3 ba_cross_nor;
    vec3_mul_cross(ba_cross_nor, ba, nor);
    vec3 cb_cross_nor;
    vec3_mul_cross(cb_cross_nor, cb, nor);
    vec3 ac_cross_nor;
    vec3_mul_cross(ac_cross_nor, ac, nor);
    vec3 ba_scaled;
    vec3_scale(ba_scaled, ba, clamp(vec3_mul_inner(ba,pa)/dot2(ba),0.0,1.0));
    vec3 cb_scaled;
    vec3_scale(cb_scaled, cb, clamp(vec3_mul_inner(cb,pb)/dot2(cb),0.0,1.0));
    vec3 ac_scaled;
    vec3_scale(ac_scaled, ac, clamp(vec3_mul_inner(ac,pc)/dot2(ac),0.0,1.0));
    vec3 ba_scaled_minus_pa;
    vec3_sub(ba_scaled_minus_pa, ba_scaled, pa);
    vec3 cb_scaled_minus_pb;
    vec3_sub(cb_scaled_minus_pb, cb_scaled, pb);
    vec3 ac_scaled_minus_pc;
    vec3_sub(ac_scaled_minus_pc, ac_scaled, pc);

    float ba_dot = vec3_mul_inner(ba_cross_nor, pa);
    float cb_dot = vec3_mul_inner(cb_cross_nor, pb);
    float ac_dot = vec3_mul_inner(ac_cross_nor, pc);

    float sum_signs = sign(ba_dot) + sign(cb_dot) + sign(ac_dot);
    bool condition = (sum_signs < 2.0);

    float ba_len2 = dot2(ba_scaled_minus_pa);
    float cb_len2 = dot2(cb_scaled_minus_pb);
    float ac_len2 = dot2(ac_scaled_minus_pc);

    float min_len2 = fmin(fmin(ba_len2, cb_len2), ac_len2);

    float nor_dot_pa = vec3_mul_inner(nor, pa);
    float proj_len2 = (nor_dot_pa * nor_dot_pa) / dot2(nor);

    float result = sqrtf(condition ? min_len2 : proj_len2);
    return result;
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

inline float sdfPlane(const vec3 pos, const vec3 normal, float height) {
    float dist = vec3_mul_inner(pos, normal) + height;
    return dist;
}

float sdfVerticalCapsule(const vec3 origin, float radius, float height, mat4x4 transf) {
    vec3 pos;
    vec3_dup(pos, origin);

    vec4 pos_homogeneous = {pos[0], pos[1], pos[2], 1.f};
    vec4 pos_rotated_homogeneous;
    mat4x4_mul_vec4(pos_rotated_homogeneous, transf, pos_homogeneous);
    vec3 pos_rotated = {
        pos_rotated_homogeneous[0],
        pos_rotated_homogeneous[1],
        pos_rotated_homogeneous[2]
    };

    pos_rotated[1] -= clamp(pos_rotated[1], 0.0f, height);

    return vec3_len(pos_rotated) - radius;
}

static float sdfCylinder(const vec3 pos,float radius, float height, mat4x4 transf) {
    vec4 pos_homogeneous = {pos[0], pos[1], pos[2], 1.f};
    vec4 pos_rotated_homogeneous;
    mat4x4_mul_vec4(pos_rotated_homogeneous, transf, pos_homogeneous);
    vec3 pos_rotated = {
        pos_rotated_homogeneous[0],
        pos_rotated_homogeneous[1],
        pos_rotated_homogeneous[2]
    };

    vec2 xz;
    xz[0] = pos_rotated[0];
    xz[1] = pos_rotated[2];
    float xzLen = vec2_len(xz);

    vec2 d;
    d[0] = xzLen - radius;
    d[1] = fabsf(pos_rotated[1]) - height;

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

const long basic_grid_dim = 10;
const long grid_dim[3] = {basic_grid_dim, basic_grid_dim, basic_grid_dim};
const long grid_strides[3] = {basic_grid_dim * basic_grid_dim, basic_grid_dim, 1};

struct SceneParams {
    float xyz[3];
    float grid[basic_grid_dim][basic_grid_dim][basic_grid_dim];
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

void scene_params_init(SceneParams *params, const SceneContext *ctx) {
    (void)ctx;
    float *raw_params = float_pointer_from_params(params);
    for (int i = 0; i < number_of_scene_params; i++) {
        raw_params[i] = 0.f;
    }
    sample_existing_sdf(grid_strides, grid_dim, &params->grid[0][0][0], [&](const vec3 pos) {
        // float ret = sdfVerticalCapsule(pos, 0.8f, 0.6f);
        // float ret = sdfSphere(pos, 0.8f);
        vec3 local;
        vec3 diff;
        vec3_set(diff, 0.3f, 0.0, 0.0f);
        vec3_sub(local, pos, diff);
        mat4x4 trans;
        mat4x4_identity(trans);
        float ret = sdfCylinder(local, 0.7f, 0.7f, trans);
        return ret;
    });
}

float scene_consistency_loss(const SceneParams *params) {
    (void)params;
    return grid_consistency_loss(grid_strides, grid_dim, &params->grid[0][0][0]);
}

void scene_consistency_gradient(const SceneParams *params, SceneParams *gradient_out) {
    (void)params;
    scene_params_fill(gradient_out, 0.0);
    grid_consistency_loss_diff(&gradient_out->grid[0][0][0], grid_strides, grid_dim, &params->grid[0][0][0]);
}

static inline float sdfFiniteSquare(const vec3 pos, float thickness) {
    // The square is centered at (0, 1, 0) and extends from -0.5 to 0.5 in x and y.
    // It has a thickness along the z-axis.

    float distToRight = pos[0] - 0.5f;
    float distToLeft = -pos[0] - 0.5f;
    float distToTop = pos[1] + 0.5f + 0.5f;
    float distToBottom = -pos[1] - 0.5f - 0.5f;
    float distToFront = pos[2] - thickness * 0.5f;
    float distToBack = -pos[2] - thickness * 0.5f;

    // The SDF of the finite square is the maximum of the distances to all six planes
    // (four for the square boundary in xy, and two for the thickness in z).
    return fmaxf(fmaxf(fmaxf(distToRight, distToLeft), fmaxf(distToTop, distToBottom)), fmaxf(distToFront, distToBack));
}

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
static inline float object_rightwall_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)params;
    (void)ctx;
    vec3 plane_normal = {-1.0, 0.0, 0.0};
    return sdfPlane(pos, plane_normal, 1.0f);
}
static inline void object_rightwall_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    default_scene_sample(sample);
    vec3_set(sample->ambient, 0.14f, 0.45f, 0.091f);
    vec3_set(sample->diffuse, 0.14f, 0.45f, 0.091f);
    vec3_set(sample->specular, 0.4f);
}

// Floor:
static inline float object_bottomwall_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)pos;
    (void)params;
    (void)ctx;
    vec3 plane_normal = {0.0, -1.0, 0.0};
    return sdfPlane(pos, plane_normal, 1.0f);
}
static inline void object_bottomwall_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    default_scene_sample(sample);
    vec3_set(sample->ambient, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->diffuse, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->specular, 0.4f);
}
inline float object_capsule_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)params;
    (void)ctx;
    vec3 world = {0.f, -0.2f, 0.f};
    vec3 local;
    vec3_sub(local, pos, world);
    mat4x4 transf;
    mat4x4_identity(transf);
    mat4x4_rotate_X(transf, transf, params->xyz[0] - M_PIf / 4.f);
    mat4x4_rotate_Y(transf, transf, params->xyz[1] - M_PIf / 4.f);
    mat4x4_rotate_Z(transf, transf, params->xyz[2] + M_PIf / 4.f);
    return sdfVerticalCapsule(local, 0.2f, 0.5f, transf);
}

inline void object_capsule_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    vec3_set(sample->diffuse, 3.f, 3.f, 0.1f);
    vec3_set(sample->specular, 0.9f, 0.9f, 0.9f);
    vec3_set(sample->ambient, 0.1f, 0.1f, 0.1f);
}

inline float object_eye_0_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)params;
    (void)ctx;
    vec3 world = {-0.02f, -0.3f, 0.5f};
    vec3 local;
    vec3_sub(local, pos, world);
    mat4x4 transf;
    mat4x4_identity(transf);
    return sdfVerticalCapsule(local, 0.05f, 0.0, transf);
}

inline void object_eye_0_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    vec3_set(sample->diffuse, 0.3f, 0.3f, 0.3f);
    vec3_set(sample->specular, 0.5f, 0.5f, 0.5f);
    vec3_set(sample->ambient, 0.1f, 0.1f, 0.1f);
}

inline float object_eye_1_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)params;
    (void)ctx;
    vec3 world = {0.1f, -0.3f, 0.6f};
    vec3 local;
    vec3_sub(local, pos, world);
    mat4x4 transf;
    mat4x4_identity(transf);
    return sdfVerticalCapsule(local, 0.1f, 0.0, transf);
}

inline void object_eye_1_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    vec3_set(sample->diffuse, 0.3f, 0.3f, 0.3f);
    vec3_set(sample->specular, 0.5f, 0.5f, 0.5f);
    vec3_set(sample->ambient, 0.1f, 0.1f, 0.1f);
}

inline float object_leg_0_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)params;
    (void)ctx;
    vec3 world = {0.2f, 0.65f, 0.1f};
    vec3 local;
    vec3_sub(local, pos, world);
    mat4x4 transf;
    mat4x4_identity(transf);
    mat4x4_rotate_Z(transf, transf, M_PIf / 3.f);
    return sdfCylinder(local, 0.1f, 0.07f, transf);
}

inline void object_leg_0_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    vec3_set(sample->diffuse, 0.3f, 0.3f, 0.3f);
    vec3_set(sample->specular, 0.5f, 0.5f, 0.5f);
    vec3_set(sample->ambient, 0.1f, 0.1f, 0.1f);
}

inline float object_leg_1_sd(const vec3 pos, const SceneParams *params, const SceneContext *ctx) {
    (void)params;
    (void)ctx;
    vec3 world = {-0.14f, 0.55f, 0.2f};
    vec3 local;
    vec3_sub(local, pos, world);
    mat4x4 transf;
    mat4x4_identity(transf);
    mat4x4_rotate_X(transf, transf, M_PIf / 5.f);
    return sdfCylinder(local, 0.2f, 0.07f, transf);
}

inline void object_leg_1_mat(const vec3 pos, const SceneParams *params, const SceneContext *ctx, SdfResult *sample) {
    (void)pos;
    (void)params;
    (void)ctx;
    vec3_set(sample->diffuse, 0.3f, 0.3f, 0.3f);
    vec3_set(sample->specular, 0.5f, 0.5f, 0.5f);
    vec3_set(sample->ambient, 0.1f, 0.1f, 0.1f);
}

extern float __enzyme_autodiff_normal(void *, int, const float *, float *, int, const SceneParams *, int, const SceneContext *ctx);

#define SCENE\
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
        enzyme_const, params, enzyme_const, ctx))\
    OBJ(object_capsule_sd,   object_capsule_mat,   __enzyme_autodiff_normal(\
        (void*)object_capsule_sd,   enzyme_dup, pos, normal,\
        enzyme_const, params, enzyme_const, ctx))\
    OBJ(object_eye_0_sd,      object_eye_0_mat,      __enzyme_autodiff_normal(\
        (void*)object_eye_0_sd,      enzyme_dup, pos, normal,\
        enzyme_const, params, enzyme_const, ctx))\
    OBJ(object_eye_1_sd,      object_eye_1_mat,      __enzyme_autodiff_normal(\
        (void*)object_eye_1_sd,      enzyme_dup, pos, normal,\
        enzyme_const, params, enzyme_const, ctx))\
    OBJ(object_leg_0_sd,      object_leg_0_mat,      __enzyme_autodiff_normal(\
        (void*)object_leg_0_sd,      enzyme_dup, pos, normal,\
        enzyme_const, params, enzyme_const, ctx))\
    OBJ(object_leg_1_sd,      object_leg_1_mat,      __enzyme_autodiff_normal(\
        (void*)object_leg_1_sd,      enzyme_dup, pos, normal,\
        enzyme_const, params, enzyme_const, ctx))

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
