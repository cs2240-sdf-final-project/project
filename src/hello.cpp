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

typedef struct {
    float distance;
    vec3 normal;
} SdfResult;


void vec2_abs(vec2 out, const vec2 in) {
    for (int i = 0; i < 2; i++) {
        out[i] = fabsf(in[i]);
    }
}
void sdfCylinder(const vec3 origin,float radius, float height, float param, SdfResult *result){
    vec3 pos;
    vec3_dup(pos, origin);
    pos[2] += param;

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


    result->distance = dist;
    vec3_norm(result->normal, pos);
}
void sdfPlane(const vec3 origin, const vec3 normal, const float height, float param, SdfResult *result) {
   
    
    vec3 pos;
    vec3_dup(pos, origin);
    pos[2] += param;

    float dist = vec3_mul_inner(pos, normal) + height;
    result->distance = dist;
    vec3_norm(result->normal, pos);
}


void sdfTriPrism(const vec3 orgin, const float h0, const float h1 , float param, SdfResult *result) {
    //h[0] represents half the length of the base of the triangular prism
    //h[1] represents half the height of the prism along the z-axis
    vec3 pos;
    vec3_dup(pos, orgin);
    pos[2] += param;

    vec3 q = {fabsf(pos[0]), fabsf(pos[1]), fabsf(pos[2])};
    float dist = fmaxf(q[2] - h1, fmaxf(q[0] * 0.866025f + pos[1] * 0.5f, -pos[1]) - h0 * 0.5f);
    result->distance = dist;
    vec3_norm(result->normal, pos);
}

void sdfVerticalCapsule(const vec3 origin, float height, float radius, float param, SdfResult *result) {
    vec3 pos;
    vec3_dup(pos, origin);
    pos[2] += param;

    pos[1] -= clamp(pos[1], 0.0f, height);

    result->distance = vec3_len(pos) - radius;
    vec3_norm(result->normal, pos);
}



void sdf_sphere(const vec3 pos, float param, SdfResult *result) {
    vec3 origin;
    vec3_set(origin, 0.0);

    origin[2] += param;

    vec3 displacement;
    vec3_sub(displacement, pos, origin);

    result->distance = vec3_len(displacement) - 3.5f;
    vec3_norm(result->normal, displacement);
}

void sdf(const vec3 pos, float param, SdfResult *result) {
    //sdfCylinder(pos, 0.5, 1.0, param, result);
    //vec3 normal = {1.0, 0.0, 0.0};
    //sdfPlane(pos, normal, 0.3, 0.5, result);
    sdfVerticalCapsule(pos, 2.0, 1.0, param, result);
    //sdfTriPrism(pos, 1.0, 5.0, param, result);
    
}


float sdf_normal_wrapper(const vec3 pos, float param) {
    SdfResult res;
    sdf(pos, param, &res);
    return res.distance;
}

extern void __enzyme_autodiff_normal(void *, int, const float *, float *, int, float);

void sdf2(const vec3 pos, float param, SdfResult *result) {
    vec3 dpos;
    vec3_set(dpos, 0.0f);

    __enzyme_autodiff_normal(
        (void*)sdf_normal_wrapper,
        enzyme_dup, pos, dpos,
        enzyme_const, param
    );
    result->distance = sdf_normal_wrapper(pos, param);
    vec3_dup(result->normal, dpos);
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

            // lightDir = -normalize(lightDirections[i].xyz);
            // vec4 lightColor = lightColors[i];

            // // Diffuse term
            // fragColor += k_d * clamp(max(dot(lightDir, normal), 0.0),0.0,1.0) * cDiffuse * lightColor;
void phongLight(vec3 radiance, const vec3 normal) {

    // vec3 light_dir = {-3.f, 0.f, -2.f };
    // //vec3 light_dir = {0.f, 1.f, -0.f };
    // vec3_norm(light_dir, light_dir);

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

    // diffuse color of the object
    vec3 cDiffuse = {0.5f, 0.3f, 0.9f};

    vec3_set(radiance, 0.f);
    for (int l = 0; l < 3; l++) {
        vec3 lightColor;
        vec3_dup(lightColor, lightColors[l]);
        vec3 light_dir;
        vec3_norm(light_dir, lightDirections[l]);

        float diffuse = clamp(fmaxf(vec3_mul_inner(light_dir, normal), 0.0f), 0.0f, 1.0f);
        radiance[0] += kd * diffuse * cDiffuse[0] * lightColor[0];
        radiance[1] += kd * diffuse * cDiffuse[1] * lightColor[1];
        radiance[2] += kd * diffuse * cDiffuse[2] * lightColor[2];
    }
}

void render_get_radiance(vec3 radiance, RandomState *rng, const vec3 origin, const vec3 direction, float param) {
    vec3_set(radiance, 0.0);
    vec3 current_position;
    vec3_dup(current_position, origin);
    for (int i = 0; i < 100; i++) {
        SdfResult res;
        sdf2(current_position, param, &res);
        if (res.distance < 1e-4f) {
            vec3 color;
            //color_normal(color, res.normal);

            phongLight(color, res.normal);
            vec3_dup(radiance, color);
            break;
        } else {
            ray_step(current_position, direction, res.distance);
        }
    }
}

extern void __enzyme_fwddiff_radiance(void *, int, float *, float *, int, RandomState *, int, const vec3, int, const vec3, int, float, float);

void render_get_gradient_helper(vec3 real, vec3 gradient, RandomState *rng, const vec3 origin, const vec3 direction) {
    float param = 0.1f;
    float d_param = 1.0f;

    vec3 radiance;
    vec3_set(radiance, 1.f);
    vec3 d_radiance;
    vec3_set(d_radiance, 1.f);

    __enzyme_fwddiff_radiance(
        (void*)render_get_radiance,
        enzyme_dup, radiance, d_radiance,
        enzyme_const, rng,
        enzyme_const, origin,
        enzyme_const, direction,
        enzyme_dup, param, d_param);

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

