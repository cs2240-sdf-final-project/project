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

float directional_derivative_inner(const vec3 pos, const vec3 direction, float t, const SceneParams *params) {
    vec3 scaled;
    vec3_scale(scaled, direction, t);
    vec3 added;
    vec3_add(added, pos, scaled);

    SceneSample sample;
    scene_sample(added, params, &sample);
    return sample.distance;
}

extern float __enzyme_fwddiff_directional(void *,
    int, const float *,
    int, const float *,
    int, float, float,
    int, const SceneParams *);

float directional_derivative(const vec3 pos, const vec3 direction, const SceneParams *params) {
    float t = 0.0f;
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
    SceneSample sample;
    scene_sample(pos, &scene_params, &sample);
    return sample.distance;
}

extern void __enzyme_autodiff_normal(void *, int, const float *, float *, int, const float *);

void get_normal_from(const vec3 pos, const SceneParams *params, vec3 normal) {
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

void ray_step(vec3 origin, const vec3 direction, float t) {
    vec3 step_size;
    vec3_dup(step_size, direction);
    vec3_scale(step_size, step_size, t);
    vec3_add(origin, origin, step_size);
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

void render_get_radiance(vec3 radiance, RandomState *rng, const vec3 origin, const vec3 direction, const SceneParams *params) {
    vec3_set(radiance, 0.0);
    vec3 current_position;
    vec3_dup(current_position, origin);
    for (int i = 0; i < 100; i++) {
        SceneSample res;
        scene_sample(current_position, params, &res);
        // float directionalDeriv = directional_derivative(current_position, direction, params)

        if (res.distance < 1e-4f) {
            vec3 normal;
            get_normal_from(current_position, params, normal);
            vec3 color;
            phongLight(color, direction, normal, &res);
            vec3_dup(radiance, color);
            break;
        } else {
            ray_step(current_position, direction, res.distance);
        }
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
        char value = (char)clamp(radiance[p] * 255.f, 0.f, 255.f);
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
