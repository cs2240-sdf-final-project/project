#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstdint>

#include "util.h"
#include "params.h"

typedef struct {
    long row_stride;
    long col_stride;
    long subpixel_stride;
    long parameter_stride;
} GradientStrides;

struct GradientImage {
    GradientStrides strides;
    long image_width;
    long image_height;
    long num_subpixels;
    float *buf;
};

typedef struct {
    long row_stride;
    long col_stride;
    long subpixel_stride;
} Strides;

struct Image {
    Strides strides;
    long image_width;
    long image_height;
    long num_floats;
    float *buf;
};

inline long get_index(const Strides *s, long r, long c, long p) {
    return r * s->row_stride + c * s->col_stride + p * s->subpixel_stride;
}

void image_fetch_dims(const Image *image, long *width, long *height) {
    *width = image->image_width;
    *height = image->image_height;
}

GradientImage *make_gradient_image(long image_width, long image_height) {
    const long num_subpixels = 3;
    long num_floats = number_of_scene_params * image_width * image_height * num_subpixels;
    assert(num_floats > 0);
    float *buf = new float[(size_t)num_floats];
    assert(buf);
    for (long i = 0; i < num_floats; i++) {
        buf[i] = 0.f;
    }
    GradientStrides strides;
    strides.parameter_stride = 1;
    strides.subpixel_stride = number_of_scene_params;
    strides.col_stride = number_of_scene_params * num_subpixels;
    strides.row_stride = number_of_scene_params * num_subpixels * image_width;
    GradientImage *ret = new GradientImage;
    ret->strides = strides;
    ret->image_height = image_height;
    ret->image_width = image_width;
    ret->num_subpixels = 3;
    ret->buf = buf;
    return ret;
}

void free_gradient_image(GradientImage *image) {
    delete[] image->buf;
    delete image;
}

inline long gradient_image_get_index(const GradientStrides *s, long r, long c, long subpixel, long param) {
    return r * s->row_stride + c * s->col_stride + subpixel * s->subpixel_stride + param * s->parameter_stride;
}

void gradient_image_set(const SceneParamsPerChannel *ppc, GradientImage *image, long ir, long ic) {
    for (long subpixel = 0; subpixel < 3; subpixel++) {
        const float *raw_params = ppc_get_channel(ppc, subpixel);
        for (long p = 0; p < number_of_scene_params; p++) {
            long index = gradient_image_get_index(&image->strides, ir, ic, subpixel, p);
            image->buf[index] = raw_params[p];
        }
    }
}

void gradient_image_get(SceneParamsPerChannel *ppc, const GradientImage *image, long ir, long ic) {
    for (long subpixel = 0; subpixel < 3; subpixel++) {
        float *raw_params = ppc_get_channel(ppc, subpixel);
        for (long p = 0; p < number_of_scene_params; p++) {
            long index = gradient_image_get_index(&image->strides, ir, ic, subpixel, p);
            raw_params[p] = image->buf[index];
        }
    }
}

void image_set(Image *image, long ir, long ic, const vec3 radiance) {
    assert(ir >= 0);
    assert(ir < image->image_height);
    assert(ic >= 0);
    assert(ic < image->image_width);
    for (long p = 0; p < 3; p++) {
        long index = get_index(&image->strides, ir, ic, p);
        image->buf[index] = radiance[p];
    }
}

void gradient_image_slice(Image *image, const GradientImage *gradient, long parameter_no) {
    SceneParamsPerChannel *ppc = make_ppc();
    for (long r = 0; r < image->image_height; r++) {
        for (long c = 0; c < image->image_width; c++) {
            gradient_image_get(ppc, gradient, r, c);
            vec3 radiance;
            for (int ch = 0; ch < 3; ch++) {
                float *raw_params = ppc_get_channel(ppc, ch);
                radiance[ch] = raw_params[parameter_no];
            }
            vec3 half;
            vec3_set(half, 0.5f);
            vec3_scale(radiance, radiance, 0.5f);
            vec3_add(radiance, radiance, half);
            image_set(image, r, c, radiance);
        }
    }
    free_ppc(ppc);
}

//////////////////////////////////////
// Begin PPM Parser from ChatGPT /////
void skip_whitespace_and_comments(FILE *f) {
    int c;
    while ((c = fgetc(f)) != EOF) {
        if (c == '#') {
            // Skip the comment line
            while ((c = fgetc(f)) != '\n' && c != EOF);
        } else if (c != ' ') {
            ungetc(c, f);
            break;
        }
    }
}

int image_read_bpm(Image *image, FILE *f) {
    char header[3];
    if (fscanf(f, "%2s", header) != 1 || strcmp(header, "P6") != 0) {
        fprintf(stderr, "Unsupported or invalid PPM format\n");
        return 0;
    }
    skip_whitespace_and_comments(f);
    long found_width;
    if (fscanf(f, "%ld", &found_width) != 1) return 0;
    assert(found_width == image->image_width);

    skip_whitespace_and_comments(f);
    long found_height;
    if (fscanf(f, "%ld", &found_height) != 1) return 0;
    assert(found_height == image->image_height);

    skip_whitespace_and_comments(f);
    int maxval;
    if (fscanf(f, "%d", &maxval) != 1 || maxval != 255) {
        fprintf(stderr, "Unsupported max color value (only 255 supported)\n");
        return 0;
    }
    // Skip the single whitespace character after maxval
    fgetc(f);
    uint8_t *buf = new uint8_t[(size_t)image->num_floats];
    size_t bytes_read = fread(buf, sizeof(uint8_t), (size_t)image->num_floats, f);
    assert(bytes_read == (size_t)image->num_floats);
    for (long i = 0; i < image->num_floats; i++) {
        image->buf[i] = (float)buf[i] * (1.f / 255.f);
    }
    delete[] buf;
    return 1; // success
}
// End PPM Parser from ChatGPT ///////
//////////////////////////////////////

Image *make_image(long image_width, long image_height) {
    long num_floats = image_width * image_height * 3;
    assert(num_floats > 0);
    float *buf = new float[(size_t)num_floats];
    assert(buf);
    Strides strides;
    strides.row_stride = image_width * 3;
    strides.col_stride = 3;
    strides.subpixel_stride = 1;
    Image *image = new Image;
    image->strides = strides;
    image->image_width = image_width;
    image->image_height = image_height;
    image->buf = buf;
    image->num_floats = num_floats;
    return image;
}

void free_image(Image *image) {
    delete[] image->buf;
    delete image;
}

// https://nullprogram.com/blog/2017/11/03/
void image_write_ppm(const Image *image, FILE *f) {
    assert(f);
    fprintf(f, "P6\n%ld %ld\n255\n", image->image_width, image->image_height);
    assert(image->num_floats > 0);
    uint8_t *buf = new uint8_t[(size_t)image->num_floats];
    for (long i = 0; i < image->num_floats; i++) {
        uint8_t value = (uint8_t)(clamp(roundf(image->buf[i] * 255.f), 0.f, 255.f));
        buf[i] = value;
    }
    fwrite(buf, 1, (size_t)image->num_floats, f);
    delete[] buf;
    fflush(f);
}

void image_get(vec3 radiance, const Image *image, long ir, long ic) {
    assert(ir >= 0);
    assert(ir < image->image_height);
    assert(ic >= 0);
    assert(ic < image->image_width);
    for (long p = 0; p < 3; p++) {
        long index = get_index(&image->strides, ir, ic, p);
        radiance[p] = image->buf[index];
    }
}

void image_copy(const Image *src, Image *dst) {
    for (long ir = 0; ir < src->image_height; ir++) {
        for (long ic = 0; ic < src->image_width; ic++) {
            vec3 radiance;
            image_get(radiance, src, ir, ic);
            image_set(dst, ir, ic, radiance);
        }
    }
}
