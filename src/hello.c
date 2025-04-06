// test.c

#include <stdio.h>
#include <assert.h>
#include "linmath.h"
#include "sim_random.h"

extern double __enzyme_autodiff(void*, double);

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

typedef struct {
    float distance;
    vec3 normal;
} SdfResult;

SdfResult sdf(const vec3 pos) {
    vec3 origin;
    vec3_set(origin, 0.0);

    vec3 displacement;
    vec3_sub(displacement, pos, origin);

    SdfResult result;
    result.distance = vec3_len(displacement);
    vec3_norm(result.normal, displacement);
    return result;
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

void render_get_radiance(vec3 radiance, RandomState *rng, vec3 origin, const vec3 direction) {
    vec3_set(radiance, 0.0);
    for (int i = 0; i < 100; i++) {
        SdfResult res = sdf(origin);
        if (res.distance < 0.0) {
            color_normal(radiance, res.normal);
            break;
        }
        ray_step(origin, direction, res.distance);
    }
}


typedef struct {
    long row_stride;
    long col_stride;
    long subpixel_stride;
} Strides;

typedef struct {
    long image_width;
    long image_height;
    float *image_out;
    Strides strides;
} Image;

long get_index(Strides *s, long r, long c, long p) {
    return r * s->row_stride + c * s->col_stride + p * s->subpixel_stride;
}

void image_set(Image *image, long ir, long ic, const vec3 radiance) {
    assert(ir >= 0);
    assert(ir < image->image_height);
    assert(ic >= 0);
    assert(ic < image->image_width);
    for (long p = 0; p < 3; p++) {
        long index = get_index(&image->strides, ir, ic, p);
        image->image_out[index] = radiance[p];
    }
}

void render_image(Image *image, RandomState *rng) {
    float aspect = image->image_height / image->image_width;
    float near_clip = 0.1;
    float far_clip = 100.0;
    float y_fov = 1.0;
    vec3 camera_position = {-10.0, 0.0, 0.0};
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

    for (long ir = 0; ir < image->image_height; ir++) {
        for (long ic = 0; ic < image->image_width; ic++) {
            float r = (float)ir;
            float c = (float)ic;

            float device_x = lerp(ic, 0.0, (float)image->image_width, -1.0, 1.0);
            float device_y = lerp(ir, 0.0, (float)image->image_height, -1.0, 1.0);

            vec4 unprojected = {device_x, device_y, 1.0, 1.0};
            vec4 homo;
            mat4x4_mul_vec4(homo, world_from_camera, unprojected);
            vec3 far;
            dehomogenize(far, homo);

            vec3 direction;
            vec3_sub(direction, far, camera_position);
            vec3_norm(direction, direction);

            vec3 radiance;
            render_get_radiance(radiance, rng, camera_position, direction);
            image_set(image, ir, ic, radiance);
        }
    }
}


double square(double x) {
    vec3 foo = {1.0, 1.0, 1.0};
    vec3 ret;
    vec3_scale(ret, foo, x);
    return vec3_len(ret);
}

double dsquare(double x) {
    // This returns the derivative of square or 2 * x
    return __enzyme_autodiff((void*) square, x);
}

int main() {
    for(double i=1; i<5; i++)
        printf("square(%f)=%f, dsquare(%f)=%f", i, square(i), i, dsquare(i));
    return 0;
}
