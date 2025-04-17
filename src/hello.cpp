#include <stdlib.h>
#include <stdio.h>
#include <stdio.h>
#include <assert.h>
#include "linmath.h"
#include "sim_random.h"

int enzyme_dup;
int enzyme_dupnoneed;
int enzyme_out;
int enzyme_const;

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

void vec3_clamp(vec3 out, const vec3 x, const vec3 min, const vec3 max) {
    vec3_min(out, x, max);
    vec3_max(out, out, min);
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


typedef struct {
    float offset;
} SceneParams;

void params_from_float_pointer(const float *params, SceneParams *out) {
    out->offset = params[0];
}

const float *float_pointer_from_params(const SceneParams *out) {
    return &out->offset;
}

typedef struct {
    float distance;
    float ambient[3];
    float diffuse[3];
    float specular[3];
    float shininess;
} SdfResult;

static inline void default_scene_sample(SdfResult *s) {
    s->distance = INFINITY;
    vec3_set(s->ambient, 0.0);
    vec3_set(s->diffuse, 0.0);
    vec3_set(s->specular, 0.0);
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
    destination->shininess = b->shininess;
}

// static inline void object_foreground_capsule(const vec3 pos, const SceneParams *params, SdfResult *sample) {
//     vec3 offset = {-1.0, 1.0, 6.0};
//     offset[0] += params->offset;
//     vec3 pos1;
//     vec3_add(pos1, pos, offset);
//     sample->distance = sdfCylinder(pos1, 1.0f, 2.0f);
//     vec3_set(sample->ambient, 0.1f);
//     vec3 cDiffuse = {0.3f, 0.5f, 0.8f};
//     vec3_dup(sample->diffuse, cDiffuse);
// }


// static inline void object_foreground(const vec3 pos, const SceneParams *params, SdfResult *sample) {
//     sample->distance = sdfVerticalCapsule(pos, 2.0, 1.0);
//     vec3_set(sample->diffuse, 1.0f);
//     vec3_set(sample->specular, 0.1f);
// }

// inline void scene_sample(const vec3 pos, const SceneParams *params, SdfResult *sample) {
//     default_scene_sample(sample);

//     SdfResult working;

//     default_scene_sample(&working);
//     object_foreground_capsule(pos, params, &working);
//     compose_scene_sample(sample, &working);

//     default_scene_sample(&working);
//     object_foreground(pos, params, &working);
//     compose_scene_sample(sample, &working);
// }
static inline void object_cylinder(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    vec3 offset = {-0.4f, -0.2f, 0.2f};
    offset[0] += params->offset;
    vec3 pos_offset;
    vec3_add(pos_offset, pos, offset);
    sample->distance = sdfCylinder(pos_offset, 0.3f, 0.9f);
    vec3_set(sample->diffuse, 0.4860f, 0.6310f, 0.6630f);
    vec3_set(sample->ambient, 0.4860f, 0.6310f, 0.6630f);
    vec3_set(sample->specular, 0.8f, 0.8f, 0.8f);
}
static inline void object_capsule(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    vec3 offset = {0.4f, -0.3f, -0.5f};
    vec3 pos_offset;
    vec3_add(pos_offset, pos, offset);
    sample->distance = sdfVerticalCapsule(pos_offset, 0.5f, 0.3f);
    vec3_set(sample->diffuse, 0.4860f, 0.6310f, 0.6630f);
    vec3_set(sample->ambient, 0.4860f, 0.6310f, 0.6630f);
    vec3_set(sample->specular, 0.8f, 0.8f, 0.8f);
}
// Back:
static inline void object_backwall(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    vec3 plane_normal = {0.0, 0.0, 1.0};
    sample->distance = sdfPlane(pos, plane_normal, 1.0f);
    vec3_set(sample->ambient, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->diffuse, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->specular, 0.4f);
}
// Ceiling:
static inline void object_topwall(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    vec3 plane_normal = {0.0, 1.0, 0.0};
    sample->distance = sdfPlane(pos, plane_normal, 1.0f);
    vec3_set(sample->ambient, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->diffuse, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->specular, 0.4f);
}
// Left:
static inline void object_leftwall(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    vec3 plane_normal = {1.0, 0.0, 0.0};
    sample->distance = sdfPlane(pos, plane_normal, 1.0f);
    vec3_set(sample->ambient, 0.63f, 0.065f, 0.05f);
    vec3_set(sample->diffuse, 0.63f, 0.065f, 0.05f);
    vec3_set(sample->specular, 0.4f);
}
// Right:
static inline void object_rightwall(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    vec3 plane_normal = {-1.0, 0.0, 0.0};
    sample->distance = sdfPlane(pos, plane_normal, 1.0f);
    vec3_set(sample->ambient, 0.14f, 0.45f, 0.091f);
    vec3_set(sample->diffuse, 0.14f, 0.45f, 0.091f);
    vec3_set(sample->specular, 0.4f);
}
// Floor:
static inline void object_bottomwall(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    vec3 plane_normal = {0.0, -1.0, 0.0};
    sample->distance = sdfPlane(pos, plane_normal, 1.0f);
    vec3_set(sample->ambient, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->diffuse, 0.725f, 0.71f, 0.68f);
    vec3_set(sample->specular, 0.4f);
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
    // width of the scene band
    const float distance_threshold = 2e-1f;
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
    SceneParams scene_params;
    params_from_float_pointer(params, &scene_params);
    return scene_sample_sdf(pos, &scene_params);
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

void phongLight(vec3 radiance, const vec3 looking, const vec3 normal, const SdfResult *sample) {
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
        vec3 light_dir;
        vec3_norm(light_dir, lightDirections[l]);
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

void get_critical_point_along(
    vec3 radiance,
    const vec3 origin,
    const vec3 direction,
    const SceneParams *params,
    RandomState *rng
) {
    // we take steps at most this size in order to avoid missing
    // sign changes in the directional derivatives
    const float max_step_size = 1.0f;
    // take this many steps to balance performance with exploring the entire scene
    const int number_of_steps = 1'000;
    // if our sdf gets smaller than this amount, we will consider it an intersection with the surface
    const float contact_threshold = 1e-4f;
    float t = 0.0;
    float previous_t = 0.0;
    SearchResult critical_point = { false, 0.0 };
    IntersectionResult intersection = { false, 0.0 };
    for (int i = 0; i < number_of_steps; i++) {
        float dir_previous_t = directional_derivative(origin, direction, previous_t, params);
        float dir_t = directional_derivative(origin, direction, t, params);
        // look for a sign change in the directional derivative
        if (!critical_point.found_critical_point && ((dir_previous_t < 0) && (dir_t > 0))) {
            // let's try to find critical_point between t and previous_t
            SearchResult local_res;
            search_for_critical_point(origin, direction, params, previous_t, t, &local_res);
            if (local_res.found_critical_point) {
                critical_point = local_res;
            }
        }
        vec3 pos;
        ray_step(pos, origin, direction, t);
        float distance = scene_sample_sdf(pos, params);
        if (distance < contact_threshold) {
            intersection.found_intersection = true;
            intersection.intersection_t = t;
            break;
        }
        float step_size = fminf(max_step_size, distance);
        previous_t = t;
        t += step_size;
    }
    if (intersection.found_intersection) {
        vec3 current_position;
        ray_step(current_position, origin, direction, intersection.intersection_t);
        SdfResult sample;
        scene_sample(current_position, params, &sample);
        vec3 normal;
        get_normal_from(normal, current_position, params);
        phongLight(radiance, direction, normal, &sample);
    } else {
        // uncomment for debug colors
        // vec3 red = {1.0f, 0.0f, 0.0f};
        vec3 black = {0.0f, 0.0f, 0.0f};
        float *which;
        if (critical_point.found_critical_point) {
            which = black; // make this red for debug colors
        } else {
            which = black;
        }
        vec3_dup(radiance, which);
    }
}

void render_get_radiance(
    vec3 radiance,
    const vec3 origin,
    const vec3 direction,
    const SceneParams *params,
    RandomState *rng
) {
    vec3_set(radiance, 0.0);
    vec3 current_position;
    vec3_dup(current_position, origin);
    for (int i = 0; i < 100; i++) {
        float distance = scene_sample_sdf(current_position, params);
        if (distance < 1e-4f) {
            SdfResult sample;
            scene_sample(current_position, params, &sample);
            vec3 normal;
            get_normal_from(normal, current_position, params);
            phongLight(radiance, direction, normal, &sample);
            break;
        } else {
            ray_step(current_position, current_position, direction, distance);
        }
    }
}

void render_get_radiance_wrapper(vec3 radiance, RandomState *rng, const vec3 origin, const vec3 direction, const float *raw_params) {
    SceneParams params;
    params_from_float_pointer(raw_params, &params);
    get_critical_point_along(radiance, origin, direction, &params, rng);
}

extern void __enzyme_fwddiff_radiance(void *, int, float *, float *, int, RandomState *, int, const vec3, int, const vec3, int, const float *, const float *);

void render_get_gradient_helper(vec3 real, vec3 gradient, RandomState *rng, const vec3 origin, const vec3 direction) {
    SceneParams params;
    params.offset = 0.1f;
    const float *raw_params = float_pointer_from_params(&params);
    float d_param = 1.0f;

    vec3 radiance;
    vec3_set(radiance, 1.f);
    vec3 d_radiance;
    vec3_set(d_radiance, 1.f);

    __enzyme_fwddiff_radiance(
        (void*)render_get_radiance_wrapper,
        enzyme_dup, radiance, d_radiance,
        enzyme_const, rng,
        enzyme_const, origin,
        enzyme_const, direction,
        enzyme_dupnoneed, raw_params, &d_param);

    vec3_dup(real, radiance);
    vec3_dup(gradient, d_radiance);
}

typedef struct {
    long row_stride;
    long col_stride;
    long subpixel_stride;
} Strides;

long get_index(Strides *s, long r, long c, long p) {
    return r * s->row_stride + c * s->col_stride + p * s->subpixel_stride;
}

typedef struct {
    Strides strides;
    long image_width;
    long image_height;
    long num_bytes;
    char *buf;
} Image;

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
void image_write_bpm(Image *image, FILE *f) {
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

void render_image(Image *real, Image *gradient, RandomState *rng) {
    float aspect = (float)real->image_width / (float)real->image_height ;
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

    for (long ir = 0; ir < real->image_height; ir++) {
        if (ir % 50 == 0) {
            printf("on row %ld\n", ir);
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

            vec3 out_real;
            vec3 out_gradient;
            SceneParams params;
            params.offset = 0.0f;
            render_get_gradient_helper(out_real, out_gradient, rng, camera_position, direction);

            image_set(real, ir, ic, out_real);
            vec3_scale(out_gradient, out_gradient, 0.5);
            vec3 half = {0.5, 0.5, 0.5};
            vec3_add(out_gradient, out_gradient, half);
            image_set(gradient, ir, ic, out_gradient);
        }
    }
}

int main(int argc, char *argv[]) {
    FILE *freal = fopen("real.bpm", "w");
    FILE *fgradient = fopen("gradient.bpm", "w");

    long image_width = 500;
    long image_height = 500;
    Image real = make_image(image_width, image_height);
    Image gradient = make_image(image_width, image_height);

    RandomState rng = make_random();
    render_image(&real, &gradient, &rng);
    image_write_bpm(&real, freal);
    image_write_bpm(&gradient, fgradient);
    free_image(&real);
    free_image(&gradient);
}
