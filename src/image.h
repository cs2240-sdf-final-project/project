#ifndef IMAGE_H
#define IMAGE_H

#include "params.h"
#include "stdio.h"
#include "linmath.h"

typedef struct GradientImage GradientImage;
typedef struct Image Image;

Image *make_image(long image_width, long image_height);
void free_image(Image *image);
int image_read_bpm(Image *image, FILE *f);
void image_write_ppm(const Image *image, FILE *f);
void image_set(Image *image, long ir, long ic, const vec3 radiance);
void image_get(vec3 radiance, const Image *image, long ir, long ic);
void image_fetch_dims(const Image *image, long *width, long *height);
void image_copy(const Image *src, Image *dst);

GradientImage *make_gradient_image(long image_width, long image_height);
void free_gradient_image(GradientImage *image);
void gradient_image_set(const SceneParamsPerChannel *ppc, GradientImage *image, long ir, long ic);
void gradient_image_get(SceneParamsPerChannel *ppc, const GradientImage *image, long ir, long ic);
void gradient_image_slice(Image *image, const GradientImage *gradient, long parameter_no);

#endif // IMAGE_H
