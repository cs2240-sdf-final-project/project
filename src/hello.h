#ifndef HELLO_H
#define HELLO_H

#include <stdio.h>
#include "linmath.h"
#include "sim_random.h"

typedef struct {
    float offset;
} SceneParams;

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
void scene_params_elementwise_add(SceneParams *out_params, const SceneParams *a, const SceneParams *b);
void scene_params_elementwise_mul(SceneParams *out_params, const SceneParams *a, const SceneParams *b);
void scene_params_fill(SceneParams *params, float fill_with);

void params_increment(SceneParams &params);
Image make_image(long image_width, long image_height);
void free_image(Image *image);
int image_read_bpm(Image *image, FILE *f);
void image_write_bpm(Image *image, FILE *f);
void image_set(Image *image, long ir, long ic, const vec3 radiance);
void image_get(vec3 radiance, Image *image, long ir, long ic);
void render_image(Image *real, Image *gradient, RandomState *rng, const SceneParams *params);

typedef struct {
    SceneParams *r, *g, *b;
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
void gradient_image_get(SceneParamsPerChannel *ppc, GradientImage *image, long ir, long ic);

#endif // HELLO_H
