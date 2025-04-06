// test.c

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

void sdf(const vec3 pos, float *param, SdfResult *result) {
    vec3 origin;
    vec3_set(origin, 0.0);

    origin[1] += *param;

    vec3 displacement;
    vec3_sub(displacement, pos, origin);

    result->distance = vec3_len(displacement) - 3.5;
    vec3_norm(result->normal, displacement);
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

void render_get_radiance(vec3 radiance, RandomState *rng, const vec3 origin, const vec3 direction, float *param) {
    vec3 current_position;
    vec3_dup(current_position, origin);
    vec3_set(radiance, 0.0);
    float extension = 3e-1;
    for (int i = 0; i < 100; i++) {
        SdfResult res;
        sdf(current_position, param, &res);
        if (res.distance > extension) {
            ray_step(current_position, direction, res.distance);
        } else {
            float intersect_prob = lerp(res.distance, -extension, extension, 1.0, 0.0);
            int should_intersect = random_next_float(rng) < intersect_prob;
            if (should_intersect) {
                vec3 color;
                color_normal(color, res.normal);
                vec3_scale(color, color, 1.f / intersect_prob);
                vec3_dup(radiance, color);
                break;
            } else {
                break;
            }
        }
    }
}

float render_get_radiance_helper(RandomState *rng, const vec3 origin, const vec3 direction, float *param, int which) {
    vec3 out;
    render_get_radiance(out, rng, origin, direction, param);
    return out[which];
}

void render_get_radiance_helper_wrapped(float *out, RandomState *rng, const vec3 origin, const vec3 direction, float *param, int which) {
    *out = render_get_radiance_helper(rng, origin, direction, param, which);
}

extern void __enzyme_autodiff(void *, ...);

float render_get_gradient_helper(RandomState *rng, const vec3 origin, const vec3 direction, int which) {
    float param = 0.1;
    float d_param = 1.0;
/*    return render_get_radiance_helper(rng, origin, direction, &param, which);*/
/**/
/*    __enzyme_autodiff(render_get_radiance_helper,*/
/*        enzyme_const, rng,*/
/*        enzyme_const, origin,*/
/*        enzyme_const, direction,*/
/*        enzyme_dup, &param, &d_param,*/
/*        enzyme_const, which);*/


    float out = 0.0;
    float d_out = 1.0;

    __enzyme_autodiff((void*)render_get_radiance_helper_wrapped,
        enzyme_dupnoneed, &out, &d_out,
        enzyme_const,     rng,
        enzyme_const,     origin,
        enzyme_const,     direction,
        enzyme_dup,       &param, &d_param,
        enzyme_const,     which);

    return param;
}

void render_get_gradient_wrapper(vec3 out, RandomState *rng, const vec3 origin, const vec3 direction) {
    int count = 2;
    for (int i = 0; i < count; i++) {
        vec3 radiance;
        radiance[0] = render_get_gradient_helper(rng, origin, direction, 0);
        radiance[1] = render_get_gradient_helper(rng, origin, direction, 1);
        radiance[2] = render_get_gradient_helper(rng, origin, direction, 2);
        vec3_add(out, out, radiance);
    }
    vec3_scale(out, out, 1.f / (float)count);
}

typedef struct {
    long row_stride;
    long col_stride;
    long subpixel_stride;
} Strides;

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
    char *buf = malloc(num_bytes);
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

void image_write_bpm(Image *image, FILE *f) {
    fprintf(f, "P6\n%ld %ld\n255\n", image->image_width, image->image_height);
    fwrite(image->buf, 1, image->num_bytes, f);
    fflush(f);
}

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
        char value = (char)(radiance[p] * 255.0);
        image->buf[index] = value;
    }
}

void render_image(Image *image, RandomState *rng) {
    float aspect = (float)image->image_width / (float)image->image_height ;
    float near_clip = 0.1f;
    float far_clip = 100.0f;
    float y_fov = 1.0f;
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
        printf("on row %ld\n", ir);
        for (long ic = 0; ic < image->image_width; ic++) {
            float r = (float)ir;
            float c = (float)ic;

            float device_x = lerp(c, 0.0, (float)image->image_width, -1.0, 1.0);
            float device_y = lerp(r, 0.0, (float)image->image_height, -1.0, 1.0);

            vec4 unprojected = {device_x, device_y, 1.0, 1.0};
            vec4 homo;
            mat4x4_mul_vec4(homo, world_from_camera, unprojected);
            vec3 far;
            dehomogenize(far, homo);

            vec3 direction;
            vec3_sub(direction, far, camera_position);
            vec3_norm(direction, direction);

            vec3 radiance;
            render_get_gradient_wrapper(radiance, rng, camera_position, direction);
            image_set(image, ir, ic, radiance);
        }
    }
}

struct double2 {
    double x, y;
};

void f(float *out, double x) { 
    out[0] = 2.0 * x;
    out[1] = 3.0 * x;
    out[2] = 4.0 * x;
}

extern void __enzyme_fwddiff(void *, ...);
extern struct double2 __enzyme_autodiff_double(void *, ...);

int main(int argc, char *argv[]) {

    vec3 outs = { 0.0, 0.0, 0.0 };
    vec3 d_outs = { 1.0, 1.0, 1.0 };
    double x = 5.0;
    double dx = 1.0;
    
    __enzyme_fwddiff((void*)f, enzyme_dup, &outs, &d_outs, enzyme_dup, x, dx); 
    printf("%g %g %g", d_outs[0], d_outs[1], d_outs[2]);

/*    long image_width = 500;*/
/*    long image_height = 500;*/
/*    Image image = make_image(image_width, image_height);*/
/*    RandomState rng = make_random();*/
/*    render_image(&image, &rng);*/
/*    FILE *f = fopen("hello.bpm", "w");*/
/*    image_write_bpm(&image, f);*/
}



