#ifndef HELLO_H
#define HELLO_H

#include <stdio.h>
#include "linmath.h"
#include "sim_random.h"

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

Image make_image(long image_width, long image_height);
void free_image(Image *image);
int image_read_bpm(Image *image, FILE *f);
void image_write_bpm(Image *image, FILE *f);
void image_set(Image *image, long ir, long ic, const vec3 radiance);
void image_get(vec3 radiance, Image *image, long ir, long ic);
void render_image(Image *real, Image *gradient, RandomState *rng, const SceneParams *params);

#endif // HELLO_H
