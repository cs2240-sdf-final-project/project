#include <stdlib.h>
#include <stdio.h>
#include <stdio.h>
#include <assert.h>
#include "linmath.h"
#include "scene.c"
#include "sim_random.h"

int enzyme_dup;
int enzyme_dupnoneed;
int enzyme_out;
int enzyme_const;

float lerp(float x, float in_min, float in_max, float out_min, float out_max) {
    return out_min + (x - in_min) * (out_max - out_min) / (in_max - in_min);
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

void ray_step(vec3 out, const vec3 origin, const vec3 direction, float t) {
    vec3 step_size;
    vec3_dup(step_size, direction);
    vec3_scale(step_size, step_size, t);
    vec3_add(out, origin, step_size);
}

void color_normal(vec3 radiance, const vec3 normal) {
    vec3_dup(radiance, normal);
    vec3_scale(radiance, radiance, 0.5);
    vec3 to_add;
    vec3_set(to_add, 0.5);
    vec3_add(radiance, radiance, to_add);
}

void phongLight(vec3 radiance, const vec3 looking, const vec3 normal, const SceneSample *sample) {
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
    if (denom < 1e-5) {
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
    const float preliminary_distance_threshold = 5e-1;
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

inline SearchResult search_for_critical_point(const vec3 origin, const vec3 direction, const SceneParams *params, float t_min, float t_max) {
    // width of the scene band
    const float distance_threshold = 1e-1f;

    // we consider directional derivatives this large to be zero
    const float directional_derivative_threshold = 1e-4f;

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
        SearchResult ret = { false, 0.0 };
        return ret;
    }

    float best_dir_t = directional_derivative(origin, direction, best_t, params);
    if (fabsf(best_dir_t) > directional_derivative_threshold) {
        SearchResult ret = { false, 0.0 };
        return ret;
    }
    SearchResult ret = { true, best_t };
    return ret;
}

typedef struct {
    bool found_intersection;
    float intersection_t;
} IntersectionResult;

float get_critial_point_along(
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
    const float contact_threshold = 1e-4;

    float t = 0.0;
    float dir_t = -1.0;
    float previous_t;
    float dir_previous_t;

    SearchResult critical_point = { false, 0.0 };
    IntersectionResult intersection;

    for (int i = 0; i < number_of_steps; i++) {
        float dir_t = directional_derivative(origin, direction, t, params);

        // look for a sign change in the directional derivative
        if (!critical_point.found_critical_point && ((dir_t < 0) && (dir_previous_t > 0))) {
            // let's try to find critical_point between t and previous_t
            SearchResult local_res = search_for_critical_point(origin, direction, params, previous_t, t);
            if (local_res.found_critical_point) {
                critical_point = local_res;
            }
        }

        float sdf = scene_sample_sdf(current_position, params);

        if (distance < contact_threshold) {
            intersection.found_intersection = true;
            intersection.intersection_t = t;
            break;
        }

        float step_size = fminf(max_step_size, distance);
        t += step_size;
        previous_t = t;
        dir_previous_t = dir_t;
    }

    if (intersection.found_intersection) {
        vec3 current_position;
        ray_step(current_position, origin, direction, intersection.intersection_t);

        SceneSample sample;
        scene_sample(current_position, params, &sample);

        vec3 normal;
        get_normal_from(normal, current_position, params);

        phongLight(radiance, direction, normal, &sample);
    } else {
        vec3_set(radiance, 0.0);
    }
}

void render_get_radiance_wrapper(vec3 radiance, RandomState *rng, const vec3 origin, const vec3 direction, const float *raw_params) {
    SceneParams params;
    params_from_float_pointer(raw_params, &params);
    render_get_radiance(radiance, rng, origin, direction, &params);
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

void render_get_gradient_wrapper(vec3 out_real, vec3 out_gradient, RandomState *rng, const vec3 origin, const vec3 direction) {
    vec3_set(out_real, 0.f);
    vec3_set(out_gradient, 0.f);
    int count = 1;
    for (int i = 0; i < count; i++) {
        vec3 real;
        vec3 gradient;
        render_get_gradient_helper(real, gradient, rng, origin, direction);
        vec3_add(out_real, out_real, real);
        vec3_add(out_gradient, out_gradient, gradient);
    }
    vec3_scale(out_real, out_real, 1.f / (float)count);
    vec3_scale(out_gradient, out_gradient, 1.f / (float)count);
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
        char value = (char)clamp(radiance[] * 255.f, 0.f, 255.f);
        image->buf[index] = value;
    }
}

void render_image(Image *real, Image *gradient, RandomState *rng) {
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

    for (long ir = 0; ir < real->image_height; ir++) {
        if (ir % 50 == 0) {
            // printf("on row %ld\n", ir);
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
            render_get_gradient_wrapper(out_real, out_gradient, rng, camera_position, direction);

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
