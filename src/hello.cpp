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

int enzyme_dup;
int enzyme_dupnoneed;
int enzyme_out;
int enzyme_const;

// width of the scene band
const float distance_threshold = 2e-1f;

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
    return &out->offset;
}

static float *float_pointer_from_params(SceneParams *out) {
    return &out->offset;
}

float scene_parameter_get(const SceneParams *params, long p) {
    const float *raw_params = float_pointer_from_params(params);
    return raw_params[p];
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

static inline void object_foreground_capsule(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    vec3 offset = {-1.0, 1.0, 6.0};
    offset[0] += params->offset;
    vec3 pos1;
    vec3_add(pos1, pos, offset);
    sample->distance = sdfCylinder(pos1, 1.0f, 2.0f);
    vec3_set(sample->ambient, 0.1f);
    vec3 cDiffuse = {0.3f, 0.5f, 0.8f};
    vec3_dup(sample->diffuse, cDiffuse);
}

static inline void object_foreground(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    (void)params;
    sample->distance = sdfVerticalCapsule(pos, 2.0, 1.0);
    vec3_set(sample->diffuse, 1.0f);
    vec3_set(sample->specular, 0.1f);
}

inline void scene_sample(const vec3 pos, const SceneParams *params, SdfResult *sample) {
    default_scene_sample(sample);

    SdfResult working;

    default_scene_sample(&working);
    object_foreground_capsule(pos, params, &working);
    compose_scene_sample(sample, &working);

    default_scene_sample(&working);
    object_foreground(pos, params, &working);
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

void phongLight(vec3 radiance, const vec3 looking, const vec3 normal, const SdfResult *sample) {
    float lightColors[3][3] = {
        {1.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 1.0f},
    };
    float lightDirections[3][3] = {
        {-3.f, 0.f, -2.f},
        {-3.f, 2.f, 0.f},
        {0.f, 2.f, 3.f},
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
    SearchResult *critical_point
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
    critical_point->found_critical_point = false;
    IntersectionResult intersection = { false, 0.0 };
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
    }
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

static inline void params_per_channel_fill(SceneParamsPerChannel *ppc, float fill_with) {
    for (int i = 0; i < 3; i++) {
        scene_params_fill(ppc->rgb[i], fill_with);
    }
}

inline void gradient_image_set(const SceneParamsPerChannel *ppc, GradientImage *image, long ir, long ic) {
    for (long subpixel = 0; subpixel < 3; subpixel++) {
        const SceneParams *params = ppc->rgb[subpixel];
        const float *raw_params = float_pointer_from_params(params);
        for (long p = 0; p < number_of_scene_params; p++) {
            long index = gradient_image_get_index(&image->strides, ir, ic, subpixel, p);
            image->buf[index] = raw_params[p];
        }
    }
}

inline void gradient_image_get(SceneParamsPerChannel *ppc, const GradientImage *image, long ir, long ic) {
    for (long subpixel = 0; subpixel < 3; subpixel++) {
        SceneParams *params = ppc->rgb[subpixel];
        float *raw_params = float_pointer_from_params(params);
        for (long p = 0; p < number_of_scene_params; p++) {
            long index = gradient_image_get_index(&image->strides, ir, ic, subpixel, p);
            raw_params[p] = image->buf[index];
        }
    }
}

// TODO: figure this shit out
void render_get_radiance_wrapper(
    const vec3 origin,
    const vec3 direction,
    const float *raw_params,
    SearchResult *critical_point,
    float *red,
    float *green,
    float *blue
) {
    const SceneParams *params = params_from_float_pointer(raw_params);
    vec3 radiance;
    get_critical_point_along(radiance, origin, direction, params, critical_point);
    float *rgb[3] = {red, green, blue};
    for (int i = 0; i < 3; i++) {
        *rgb[i] = radiance[i];
    }
}

extern void __enzyme_autodiff_radiance(
    void *,
    int, const float *,
    int, const float *,
    int, const float *,
    int, SearchResult *,
    int, float *, float *,
    int, float *, float *,
    int, float *, float *
);

void render_pixel(
    vec3 real,
    const vec3 origin,
    const vec3 direction,
    SceneParamsPerChannel *params_per_channel,
    const SceneParams *params
) {
    params_per_channel_fill(params_per_channel, 1.f);

    SearchResult critical_point;

    const float *raw_params = float_pointer_from_params(params);

    float *raw_out_params[3];
    for (int c =0; c < 3; c++) {
        raw_out_params[c] = float_pointer_from_params(params_per_channel->rgb[c]);
    }

    vec3 radiance;
    // Calculate radiance with autodiff
    // __enzyme_autodiff_radiance(
    //     (void*)render_get_radiance_wrapper,
    //     enzyme_const, origin,
    //     enzyme_const, direction,
    //     enzyme_out, raw_params,
    //     enzyme_const, &critical_point,
    //     enzyme_dup, &radiance[RED], raw_out_params[RED],
    //     enzyme_dup, &radiance[GREEN], raw_out_params[GREEN],
    //     enzyme_dup, &radiance[BLUE], raw_out_params[BLUE]
    // );

    if(critical_point.found_critical_point) {
        vec3 y_star;
        ray_step(y_star, origin, direction, critical_point.t_if_found_critical_point);
        SdfResult sample;
        scene_sample(y_star, params, &sample);
        vec3 normal;
        get_normal_from(normal, y_star, params);

        vec3 y_star_radiance;
        phongLight(y_star_radiance, direction, normal, &sample);

        // SceneParams* paramOut;

        // diff_sdf(y_star, paramOut, params);

        // vec3 deltaL;
        // vec3_sub(deltaL, y_star_radiance, radiance);

        // vec3 boundary_integral;
        // vec3_set(boundary_integral, 0.0f);
        // TODO: figure this out with the matrix multiply
        // vec3_scale(boundary_integral, deltaL, - vn / distance_threshold);

        // vec3_add(d_radiance, d_radiance, boundary_integral);
    }
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

/** TODO: this output shape should be correct */
void render_image(Image *real, GradientImage *gradient, const SceneParams *params) {
    float aspect = (float)real->image_width / (float)real->image_height ;
    float near_clip = 0.1f;
    float far_clip = 100.0f;
    float y_fov = 1.0f;
    vec3 camera_position = {-10.0, 0.0, 10.0};
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

    SceneParamsPerChannel ppc;
    for (int c = 0; c < 3; c++) {
        ppc.rgb[c] = make_scene_params();
    }

    for (long ir = 0; ir < real->image_height; ir++) {
        // if (ir % 50 == 0) {
        //     printf("on row %ld\n", ir);
        // }
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
            render_pixel(out_real, camera_position, direction, &ppc, params);

            image_set(real, ir, ic, out_real);
            gradient_image_set(&ppc, gradient, ir, ic);
        }
    }

    for (int c = 0; c < 3; c++) {
        free_scene_params(ppc.rgb[c]);
    }
}
