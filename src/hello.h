#ifndef HELLO_H
#define HELLO_H

#include <stdio.h>
#include "linmath.h"
#include "sim_random.h"

extern int number_of_scene_params;


typedef struct SceneParams SceneParams;

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



SceneParams *make_scene_params();
void free_scene_params(SceneParams *params);
float scene_parameter_get(const SceneParams *params, long p);
void scene_params_set(SceneParams *params, long p, float value);
void scene_params_elementwise_add(SceneParams *out_params, const SceneParams *a, const SceneParams *b);
void scene_params_elementwise_mul(SceneParams *out_params, const SceneParams *a, const SceneParams *b);
void scene_params_scale(SceneParams *out_params, const SceneParams *a, float scale_by);
void scene_params_fill(SceneParams *params, float fill_with);


void params_increment(SceneParams &params);
Image make_image(long image_width, long image_height);
void free_image(Image *image);
int image_read_bpm(Image *image, FILE *f);
void image_write_ppm(Image *image, FILE *f);
void image_set(Image *image, long ir, long ic, const vec3 radiance);
void image_get(vec3 radiance, Image *image, long ir, long ic);

const int RED = 0;
const int GREEN = 1;
const int BLUE = 2;

typedef struct {
    SceneParams *rgb[3];
} SceneParamsPerChannel;

typedef struct {
    long row_stride;
    long col_stride;
    long subpixel_stride;
    long parameter_stride;
} GradientStrides;

typedef struct {
    GradientStrides strides;
    long image_width;
    long image_height;
    long num_subpixels;
    float *buf;
} GradientImage;

GradientImage make_gradient_image(long image_width, long image_height);
void free_gradient_image(GradientImage *image);
void gradient_image_set(const SceneParamsPerChannel *ppc, GradientImage *image, long ir, long ic);
void gradient_image_get(SceneParamsPerChannel *ppc, const GradientImage *image, long ir, long ic);

void render_image(Image *real, GradientImage *gradient, const SceneParams *params, RandomState *random);
void gradient_image_slice(Image *image, const GradientImage *gradient, long parameter_no);

#endif // HELLO_H
