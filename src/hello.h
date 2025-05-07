#ifndef HELLO_H
#define HELLO_H

#include <stdio.h>
#include "linmath.h"
#include "sim_random.h"

extern int number_of_scene_params;

typedef struct SceneContext SceneContext;
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
    long num_floats;
    float *buf;
} Image;

SceneParams *uninit_scene_params();
void scene_params_init(SceneParams *params, const SceneContext *ctx);

void free_scene_params(SceneParams *params);
float scene_parameter_get(const SceneParams *params, long p);
void scene_params_set(SceneParams *params, long p, float value);
void scene_params_elementwise_add(SceneParams *out_params, const SceneParams *a, const SceneParams *b);
void scene_params_elementwise_mul(SceneParams *out_params, const SceneParams *a, const SceneParams *b);
void scene_params_scale(SceneParams *out_params, const SceneParams *a, float scale_by);
void scene_params_fill(SceneParams *params, float fill_with);


float scene_consistency_loss(const SceneParams *params);
void scene_consistency_gradient(const SceneParams *params, SceneParams *gradient_out);

SceneContext *make_scene_context();
void free_scene_context(SceneContext *ctx);

Image make_image(long image_width, long image_height);
void free_image(Image *image);
int image_read_bpm(Image *image, FILE *f);
void image_write_ppm(const Image *image, FILE *f);
void image_set(Image *image, long ir, long ic, const vec3 radiance);
void image_get(vec3 radiance, const Image *image, long ir, long ic);

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

void render_image_phong(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *random);
void render_image_tracing(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *random);
void finite_differences(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *rng);
void render_image_effects(Image *real, GradientImage *gradient, const SceneParams *params, const SceneContext *ctx, RandomState *rng);

void gradient_image_slice(Image *image, const GradientImage *gradient, long parameter_no);

// typedef struct Renderer PixelRenderer;
// PixelRenderer *make_renderer(long image_width, long image_height);
// void free_renderer(PixelRenderer *renderer);
// void project_pixel_get_gradient(vec3 real, SceneParamsPerChannel *ppc, PixelRenderer *renderer, long ir, long ic, const SceneParams *params, const SceneContext *ctx);
// void project_pixel_get_radiance(vec3 real, PixelRenderer *renderer, long ir, long ic, const SceneParams *params, const SceneContext *ctx);

#endif // HELLO_H
